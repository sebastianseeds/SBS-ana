//intime cluster selection analysis
int cluster_intime_select(int num_hcal_clusid,double hcal_clus_atime[],double atime_sh,double hcal_clus_e[],double coin_mean,double coin_sig_fac,double coin_profile_sigma){
  
  double maxE = 0.0;
  double best_cluster_index;

  //loop through all clusters and eliminate all clusters out of time, then get cluster with max energy from those.
  for(int c = 0; c<num_hcal_clusid; c++){ //assuming num_hcal_clusid is essentially Ndata.sbs.hcal.clus.id
	
    double atime = hcal_clus_atime[c];
    double atime_diff = atime - atime_sh;
    double clus_energy = hcal_clus_e[c];
    bool passCoin = abs(atime_diff - coin_mean) < coin_sig_fac*coin_profile_sigma;

    //now check if the current cluster passed the intime check
    if( passCoin ){
      //if intime is passed, check to see if the current cluster has the highest energy on this event
      bool new_max_E = maxE < clus_energy;
      
      //if the current cluster has the highest energy, update the best cluster index and overwrite the event cluster energy maximum to compare with the next cluster on this event
      if( new_max_E ){
	maxE = clus_energy;
	best_cluster_index = c;
      }
    }

  }//end for loop over clusters

  return best_cluster_index;
}
