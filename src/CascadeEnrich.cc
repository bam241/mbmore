// Implements the CascadeEnrich class
#include "CascadeEnrich.h"
#include "behavior_functions.h"
#include "enrich_functions.h"
#include "sim_init.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <vector>
#include <boost/lexical_cast.hpp>


namespace mbmore {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  CascadeEnrich::CascadeEnrich(cyclus::Context* ctx)
    : cyclus::Facility(ctx),
  feed_recipe(""),
  max_centrifuges(),
  design_feed_assay(),
  design_product_assay(),
  design_tails_assay(),
  centrifuge_velocity(485.0),
  temp(320.0),
  height(0.5),
  diameter(0.15),
  machine_feed(15),
  max_enrich(1),
  design_feed_flow(0),
  feed_commod(""),
  product_commod(""),
  tails_commod(""),
  order_prefs(true) {}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CascadeEnrich::~CascadeEnrich() {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
std::string CascadeEnrich::str() {
  std::stringstream ss;
  ss << cyclus::Facility::str()
     << " with enrichment facility parameters:"
     << " * SWU capacity: " << SwuCapacity()
     << " * Tails assay: " << tails_assay << " * Feed assay: " << FeedAssay()
     << " * Input cyclus::Commodity: " << feed_commod
     << " * Output cyclus::Commodity: " << product_commod
     << " * Tails cyclus::Commodity: " << tails_commod;
  return ss.str();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void CascadeEnrich::Build(cyclus::Agent* parent) {
  using cyclus::Material;

  tails_assay = design_tails_assay;
  
  // Calculate ideal machine performance
  double design_delU = CalcDelU(centrifuge_velocity, height, diameter,
				Mg2kgPerSec(machine_feed), temp,
				cut, eff, M, dM, x, flow_internal);
  double design_alpha = AlphaBySwu(design_delU, Mg2kgPerSec(machine_feed),
				   cut, M);

  // Design ideal cascade based on target feed flow and product assay
  std::pair<int, int> n_stages =
    FindNStages(design_alpha, design_feed_assay, design_product_assay,
		design_tails_assay);

  // TODO DELETE THIS, STAGES ARE ALREADY INTS
  // set as internal state variables
  // int truncates but we need # of stages to assure target values,
  // so if the number is 5.1 we need 6. 
  //  n_enrich_stages = int(n_stages.first) + 1;
  //  n_strip_stages = int(n_stages.second) + 1;
  n_enrich_stages = n_stages.first;
  n_strip_stages = n_stages.second;

  std::pair<int,double> cascade_info = DesignCascade(FlowPerSec(design_feed_flow),
						     design_alpha,
						     design_delU,
						     cut, max_centrifuges,
						     n_stages);

  max_feed_inventory = FlowPerMon(cascade_info.second);
  // Number of machines times swu per machine
  SwuCapacity(cascade_info.first * FlowPerMon(design_delU));

  Facility::Build(parent);
  if (initial_feed > 0) {
    inventory.Push(
      Material::Create(
        this, initial_feed, context()->GetRecipe(feed_recipe)));
  }
  
  LOG(cyclus::LEV_DEBUG2, "EnrFac") << "CascadeEnrich "
				    << " entering the simuluation: ";
  LOG(cyclus::LEV_DEBUG2, "EnrFac") << str();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void CascadeEnrich::Tick() {

 current_swu_capacity = SwuCapacity();
 
 }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void CascadeEnrich::Tock() {
  using cyclus::toolkit::RecordTimeSeries;

  LOG(cyclus::LEV_INFO4, "EnrFac") << prototype() << " used "
                                   << intra_timestep_swu_ << " SWU";
  RecordTimeSeries<cyclus::toolkit::ENRICH_SWU>(this, intra_timestep_swu_);
  LOG(cyclus::LEV_INFO4, "EnrFac") << prototype() << " used "
                                   << intra_timestep_feed_ << " feed";
  RecordTimeSeries<cyclus::toolkit::ENRICH_FEED>(this, intra_timestep_feed_);

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
std::set<cyclus::RequestPortfolio<cyclus::Material>::Ptr>
CascadeEnrich::GetMatlRequests() {
  using cyclus::Material;
  using cyclus::RequestPortfolio;
  using cyclus::Request;

  std::set<RequestPortfolio<Material>::Ptr> ports;
  RequestPortfolio<Material>::Ptr port(new RequestPortfolio<Material>());
  Material::Ptr mat = Request_();
  double amt = mat->quantity();

  if (amt > cyclus::eps_rsrc()) {
    port->AddRequest(mat, this, feed_commod);
    ports.insert(port);
  }

  return ports;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Sort offers of input material to have higher preference for more
//  U-235 content
void CascadeEnrich::AdjustMatlPrefs(
    cyclus::PrefMap<cyclus::Material>::type& prefs) {
  using cyclus::Bid;
  using cyclus::Material;
  using cyclus::Request;

  if (order_prefs == false) {
    return;
  }

  cyclus::PrefMap<cyclus::Material>::type::iterator reqit;

  // Loop over all requests
  for (reqit = prefs.begin(); reqit != prefs.end(); ++reqit) {
    std::vector<Bid<Material>*> bids_vector;
    std::map<Bid<Material>*, double>::iterator mit;
    for (mit = reqit->second.begin(); mit != reqit->second.end(); ++mit) {
      Bid<Material>* bid = mit->first;
      bids_vector.push_back(bid);
    }

    std::sort(bids_vector.begin(), bids_vector.end(), SortBids);

    // Assign preferences to the sorted vector
    double n_bids = bids_vector.size();
    bool u235_mass = 0;

    for (int bidit = 0; bidit < bids_vector.size(); bidit++) {
      int new_pref = bidit + 1;

      // For any bids with U-235 qty=0, set pref to zero.
      if (!u235_mass) {
        cyclus::Material::Ptr mat = bids_vector[bidit]->offer();
        cyclus::toolkit::MatQuery mq(mat);
        if (mq.mass(922350000) == 0) {
          new_pref = -1;
        } else {
          u235_mass = true;
        }
      }
      (reqit->second)[bids_vector[bidit]] = new_pref;
    }  // each bid
  }    // each Material Request
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void CascadeEnrich::AcceptMatlTrades(
    const std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                                cyclus::Material::Ptr> >& responses) {
  // see
  // http://stackoverflow.com/questions/5181183/boostshared-ptr-and-inheritance
  std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                        cyclus::Material::Ptr> >::const_iterator it;
  for (it = responses.begin(); it != responses.end(); ++it) {
    AddMat_(it->second);
  }
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void CascadeEnrich::AddMat_(cyclus::Material::Ptr mat) {
  // Elements and isotopes other than U-235, U-238 are sent directly to tails
  cyclus::CompMap cm = mat->comp()->atom();
  bool extra_u = false;
  bool other_elem = false;
  for (cyclus::CompMap::const_iterator it = cm.begin(); it != cm.end(); ++it) {
    if (pyne::nucname::znum(it->first) == 92) {
      if (pyne::nucname::anum(it->first) != 235 &&
          pyne::nucname::anum(it->first) != 238 && it->second > 0) {
        extra_u = true;
      }
    } else if (it->second > 0) {
      other_elem = true;
    }
  }
  if (extra_u) {
    cyclus::Warn<cyclus::VALUE_WARNING>(
        "More than 2 isotopes of U.  "
        "Istopes other than U-235, U-238 are sent directly to tails.");
  }
  if (other_elem) {
    cyclus::Warn<cyclus::VALUE_WARNING>(
        "Non-uranium elements are "
        "sent directly to tails.");
  }

  LOG(cyclus::LEV_INFO5, "EnrFac") << prototype() << " is initially holding "
                                   << inventory.quantity() << " total.";

  try {
    inventory.Push(mat);
  } catch (cyclus::Error& e) {
    e.msg(Agent::InformErrorMsg(e.msg()));
    throw e;
  }

  LOG(cyclus::LEV_INFO5, "EnrFac")
      << prototype() << " added " << mat->quantity() << " of " << feed_commod
      << " to its inventory, which is holding " << inventory.quantity()
      << " total.";
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
std::set<cyclus::BidPortfolio<cyclus::Material>::Ptr>
CascadeEnrich::GetMatlBids(cyclus::CommodMap<cyclus::Material>::type& out_requests) {
  using cyclus::Bid;
  using cyclus::BidPortfolio;
  using cyclus::CapacityConstraint;
  using cyclus::Converter;
  using cyclus::Material;
  using cyclus::Request;
  using cyclus::toolkit::MatVec;

  std::set<BidPortfolio<Material>::Ptr> ports;

  if ((out_requests.count(tails_commod) > 0) && (tails.quantity() > 0)) {
    BidPortfolio<Material>::Ptr tails_port(new BidPortfolio<Material>());

    std::vector<Request<Material>*>& tails_requests =
        out_requests[tails_commod];
    std::vector<Request<Material>*>::iterator it;
    for (it = tails_requests.begin(); it != tails_requests.end(); ++it) {
      // offer bids for all tails material, keeping discrete quantities
      // to preserve possible variation in composition
      MatVec mats = tails.PopN(tails.count());
      tails.Push(mats);
      for (int k = 0; k < mats.size(); k++) {
        Material::Ptr m = mats[k];
        Request<Material>* req = *it;
        tails_port->AddBid(req, m, this);
      }
    }
    // overbidding (bidding on every offer)
    // add an overall capacity constraint
    CapacityConstraint<Material> tails_constraint(tails.quantity());
    tails_port->AddConstraint(tails_constraint);
    LOG(cyclus::LEV_INFO5, "EnrFac") << prototype()
                                     << " adding tails capacity constraint of "
                                     << tails.capacity();
    ports.insert(tails_port);
  }

  if ((out_requests.count(product_commod) > 0) && (inventory.quantity() > 0)) {
    BidPortfolio<Material>::Ptr commod_port(new BidPortfolio<Material>());

    std::vector<Request<Material>*>& commod_requests =
        out_requests[product_commod];
    std::vector<Request<Material>*>::iterator it;
    for (it = commod_requests.begin(); it != commod_requests.end(); ++it) {
      Request<Material>* req = *it;
      Material::Ptr mat = req->target();
      double request_enrich = cyclus::toolkit::UraniumAssay(mat);
      if (ValidReq(req->target()) &&
          ((request_enrich < max_enrich) ||
           (cyclus::AlmostEq(request_enrich, max_enrich)))) {
        Material::Ptr offer = Offer_(req->target());
        commod_port->AddBid(req, offer, this);
      }
    }

    Converter<Material>::Ptr sc(new SWUConverter(FeedAssay(), tails_assay));
    Converter<Material>::Ptr nc(new NatUConverter(FeedAssay(), tails_assay));
    CapacityConstraint<Material> swu(swu_capacity, sc);
    CapacityConstraint<Material> natu(inventory.quantity(), nc);
    commod_port->AddConstraint(swu);
    commod_port->AddConstraint(natu);

    LOG(cyclus::LEV_INFO5, "EnrFac")
        << prototype() << " adding a swu constraint of " << swu.capacity();
    LOG(cyclus::LEV_INFO5, "EnrFac")
        << prototype() << " adding a natu constraint of " << natu.capacity();
    ports.insert(commod_port);
  }
  return ports;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void CascadeEnrich::GetMatlTrades(
    const std::vector<cyclus::Trade<cyclus::Material> >& trades,
    std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                          cyclus::Material::Ptr> >& responses) {
  using cyclus::Material;
  using cyclus::Trade;

  intra_timestep_swu_ = 0;
  intra_timestep_feed_ = 0;

  std::vector<Trade<Material> >::const_iterator it;
  for (it = trades.begin(); it != trades.end(); ++it) {
    double qty = it->amt;
    std::string commod_type = it->bid->request()->commodity();
    Material::Ptr response;

    // Figure out whether material is tails or enriched,
    // if tails then make transfer of material
    if (commod_type == tails_commod) {
      LOG(cyclus::LEV_INFO5, "EnrFac")
          << prototype() << " just received an order"
          << " for " << it->amt << " of " << tails_commod;
      double pop_qty = std::min(qty, tails.quantity());
      response = tails.Pop(pop_qty, cyclus::eps_rsrc());
    } else {
      LOG(cyclus::LEV_INFO5, "EnrFac")
          << prototype() << " just received an order"
          << " for " << it->amt << " of " << product_commod;
      response = Enrich_(it->bid->offer(), qty);
    }
    responses.push_back(std::make_pair(*it, response));
  }

  if (cyclus::IsNegative(tails.quantity())) {
    std::stringstream ss;
    ss << "is being asked to provide more than its current inventory.";
    throw cyclus::ValueError(Agent::InformErrorMsg(ss.str()));
  }
  if (cyclus::IsNegative(current_swu_capacity)) {
    throw cyclus::ValueError("EnrFac " + prototype() +
                             " is being asked to provide more than" +
                             " its SWU capacity.");
  }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cyclus::Material::Ptr CascadeEnrich::Enrich_(cyclus::Material::Ptr mat,
                                          double qty) {
  using cyclus::Material;
  using cyclus::ResCast;
  using cyclus::toolkit::Assays;
  using cyclus::toolkit::UraniumAssay;
  using cyclus::toolkit::SwuRequired;
  using cyclus::toolkit::FeedQty;
  using cyclus::toolkit::TailsQty;

  // get enrichment parameters
  Assays assays(FeedAssay(), UraniumAssay(mat), tails_assay);
  double swu_req = SwuRequired(qty, assays);
  double natu_req = FeedQty(qty, assays);

  // Determine the composition of the natural uranium
  // (ie. U-235+U-238/TotalMass)
  double pop_qty = inventory.quantity();
  Material::Ptr natu_matl = inventory.Pop(pop_qty, cyclus::eps_rsrc());
  inventory.Push(natu_matl);

  cyclus::toolkit::MatQuery mq(natu_matl);
  std::set<cyclus::Nuc> nucs;
  nucs.insert(922350000);
  nucs.insert(922380000);
  double natu_frac = mq.mass_frac(nucs);
  double feed_req = natu_req / natu_frac;

  // pop amount from inventory and blob it into one material
  Material::Ptr r;
  try {
    // required so popping doesn't take out too much
    if (cyclus::AlmostEq(feed_req, inventory.quantity())) {
      r = cyclus::toolkit::Squash(inventory.PopN(inventory.count()));
    } else {
      r = inventory.Pop(feed_req, cyclus::eps_rsrc());
    }
  } catch (cyclus::Error& e) {
    NatUConverter nc(FeedAssay(), tails_assay);
    std::stringstream ss;
    ss << " tried to remove " << feed_req << " from its inventory of size "
       << inventory.quantity()
       << " and the conversion of the material into natu is "
       << nc.convert(mat);
    throw cyclus::ValueError(Agent::InformErrorMsg(ss.str()));
  }

  // "enrich" it, but pull out the composition and quantity we require from the
  // blob
  cyclus::Composition::Ptr comp = mat->comp();
  Material::Ptr response = r->ExtractComp(qty, comp);
  tails.Push(r);

  current_swu_capacity -= swu_req;

  intra_timestep_swu_ += swu_req;
  intra_timestep_feed_ += feed_req;
  RecordEnrichment_(feed_req, swu_req);

  LOG(cyclus::LEV_INFO5, "EnrFac") << prototype()
                                   << " has performed an enrichment: ";
  LOG(cyclus::LEV_INFO5, "EnrFac") << "   * Feed Qty: " << feed_req;
  LOG(cyclus::LEV_INFO5, "EnrFac") << "   * Feed Assay: "
                                   << assays.Feed() * 100;
  LOG(cyclus::LEV_INFO5, "EnrFac") << "   * Product Qty: " << qty;
  LOG(cyclus::LEV_INFO5, "EnrFac") << "   * Product Assay: "
                                   << assays.Product() * 100;
  LOG(cyclus::LEV_INFO5, "EnrFac") << "   * Tails Qty: "
                                   << TailsQty(qty, assays);
  LOG(cyclus::LEV_INFO5, "EnrFac") << "   * Tails Assay: "
                                   << assays.Tails() * 100;
  LOG(cyclus::LEV_INFO5, "EnrFac") << "   * SWU: " << swu_req;
  LOG(cyclus::LEV_INFO5, "EnrFac") << "   * Current SWU capacity: "
                                   << current_swu_capacity;

  return response;
}
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void CascadeEnrich::RecordEnrichment_(double natural_u, double swu) {
  using cyclus::Context;
  using cyclus::Agent;

  LOG(cyclus::LEV_DEBUG1, "EnrFac") << prototype()
                                    << " has enriched a material:";
  LOG(cyclus::LEV_DEBUG1, "EnrFac") << "  * Amount: " << natural_u;
  LOG(cyclus::LEV_DEBUG1, "EnrFac") << "  *    SWU: " << swu;

  Context* ctx = Agent::context();
  ctx->NewDatum("Enrichments")
      ->AddVal("ID", id())
      ->AddVal("Time", ctx->time())
      ->AddVal("Natural_Uranium", natural_u)
      ->AddVal("SWU", swu)
      ->Record();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cyclus::Material::Ptr CascadeEnrich::Request_() {
  double qty = std::max(0.0, inventory.capacity() - inventory.quantity());
  return cyclus::Material::CreateUntracked(qty,
                                           context()->GetRecipe(feed_recipe));
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cyclus::Material::Ptr CascadeEnrich::Offer_(cyclus::Material::Ptr mat) {
  cyclus::toolkit::MatQuery q(mat);
  cyclus::CompMap comp;
  comp[922350000] = q.atom_frac(922350000);
  comp[922380000] = q.atom_frac(922380000);
  return cyclus::Material::CreateUntracked(
      mat->quantity(), cyclus::Composition::CreateFromAtom(comp));
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bool CascadeEnrich::ValidReq(const cyclus::Material::Ptr mat) {
  cyclus::toolkit::MatQuery q(mat);
  double u235 = q.atom_frac(922350000);
  double u238 = q.atom_frac(922380000);
  return (u238 > 0 && u235 / (u235 + u238) > tails_assay);
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  double CascadeEnrich::FeedAssay() {
  using cyclus::Material;

  if (inventory.empty()) {
    return 0;
  }
  double pop_qty = inventory.quantity();
  cyclus::Material::Ptr fission_matl =
      inventory.Pop(pop_qty, cyclus::eps_rsrc());
  inventory.Push(fission_matl);
  return cyclus::toolkit::UraniumAssay(fission_matl);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
extern "C" cyclus::Agent* ConstructCascadeEnrich(cyclus::Context* ctx) {
  return new CascadeEnrich(ctx);
}
  
}  // namespace mbmore
