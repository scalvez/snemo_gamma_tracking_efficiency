/* snemo_gamma_tracking_studies_module.h
 * Author(s)     : Xavier Garrido <garrido@lal.in2p3.fr>
 * Creation date : 2014-12-27
 * Last modified : 2014-12-27
 *
 * Copyright (C) 2014 Xavier Garrido <garrido@lal.in2p3.fr>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 *
 * Description:
 *
 * A module that estimates the gamma tracking efficiency by comparing
 * calorimeter reconstruction sequence to the true i.e. simulated seuqence.
 *
 * History:
 *
 */

#ifndef ANALYSIS_SNEMO_GAMMA_TRACKING_EFFICIENCY_MODULE_H_
#define ANALYSIS_SNEMO_GAMMA_TRACKING_EFFICIENCY_MODULE_H_ 1

// Standard libraires:
#include <set>
#include <map>
// Third party:
// - Bayeux/geomtools:
#include <geomtools/geom_id.h>

// Data processing module abstract base class
#include <dpp/base_module.h>

#include <snemo/datamodels/calibrated_data.h>

namespace mygsl {
  class histogram_pool;
}

namespace snemo {
  namespace geometry {
    class locator_plugin;
  }
}

namespace analysis {

  class snemo_gamma_tracking_efficiency_module : public dpp::base_module
  {
  public:

    /// Setting histogram pool
    void set_histogram_pool(mygsl::histogram_pool & pool_);

    /// Grabbing histogram pool
    mygsl::histogram_pool & grab_histogram_pool();

    /// Typedef for calorimeters collection
    typedef std::set<geomtools::geom_id> calo_list_type;

    /// Typedef for gamma dictionnaries
    typedef std::map<int, calo_list_type> gamma_dict_type;

    /// Constructor
    snemo_gamma_tracking_efficiency_module(datatools::logger::priority = datatools::logger::PRIO_FATAL);

    /// Destructor
    virtual ~snemo_gamma_tracking_efficiency_module();

    /// Initialization
    virtual void initialize(const datatools::properties  & setup_,
                            datatools::service_manager   & service_manager_,
                            dpp::module_handle_dict_type & module_dict_);

    /// Reset
    virtual void reset();

    void get_new_neighbours(geomtools::geom_id gid,
                            const snemo::datamodel::calibrated_data::calorimeter_hit_collection_type & cch,
                            std::vector<geomtools::geom_id>  & ccl,
                            std::vector<geomtools::geom_id>  & a_cluster);
    ///*,                                                             std::vector<geomtools::geom_id> new_calib_neighbours*/);

    /// Data record processing
    virtual process_status process(datatools::things & data_);

  protected:

    /// Give default values to specific class members.
    void _set_defaults();

    /// Identify the calorimeter blocks clusters from the 'particle_track_data' bank
    void _pre_process_clustering(const datatools::things & data_,
                                 gamma_dict_type & gammas_);

    /// Get gammas sequence from 'simulated_data' bank
    dpp::base_module::process_status _process_simulated_gammas(const datatools::things & data_,
                                                               gamma_dict_type & gammas_);

    /// Get gammas sequence from 'particle_track_data' bank
    dpp::base_module::process_status _process_reconstructed_gammas(const datatools::things & data_,
                                                                   gamma_dict_type & gammas_);

    /// Compare simulated and reconstructed gamma track length
    dpp::base_module::process_status _compute_gamma_track_length(const datatools::things & data_);

    /// Compare 2 sequences of calorimeters
    bool _compare_sequences(const gamma_dict_type & simulated_gammas_,
                            const gamma_dict_type & reconstructed_gammas_);

    /// Compare 2 sequences of calorimeters for gamma clustering only
    bool _compare_sequences_cluster(const gamma_dict_type & simulated_gammas_,
                            const gamma_dict_type & clustered_gammas_);

  private:

    // The key fields from 'event header' bank to build the histogram key:
    std::vector<std::string> _key_fields_;

    // The histogram pool :
    mygsl::histogram_pool * _histogram_pool_;

    // Locator plugin
    const snemo::geometry::locator_plugin * _locator_plugin_;

    /// Internal structure to compute efficiency
    struct efficiency_type {
      size_t nevent; //!< Total number of event processed
      size_t ntotal; //!< Total number of gammas simulated
      size_t ngamma; //!< Number of gammas simulated for each event
      size_t ngood;  //!< Number of gammas well reconstructed
      size_t nmiss;  //!< Number of gammas that do not trigger detector
      size_t ngood_event;  //!< Number of events fully and successfully reconstructed
      size_t nevent_gammas;  //!< Number of events with at least one gamma

      size_t no_gt_ngood_event;  //!< Number of events fully and successfully reconstructed
      size_t no_gt_nevent_gammas;  //!< Number of events with at least one gamma
    };

    /// Efficiency structure
    efficiency_type _efficiency_;

    /// No GT efficiency structure
    efficiency_type _no_gt_efficiency_;

    // Macro to automate the registration of the module :
    DPP_MODULE_REGISTRATION_INTERFACE(snemo_gamma_tracking_efficiency_module);
  };

} // namespace analysis

#endif // ANALYSIS_SNEMO_GAMMA_TRACKING_EFFICIENCY_MODULE_H_

// end of snemo_gamma_tracking_efficiency_module.h
/*
** Local Variables: --
** mode: c++ --
** c-file-style: "gnu" --
** tab-width: 2 --
** End: --
*/
