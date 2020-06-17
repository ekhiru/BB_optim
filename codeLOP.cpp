int kendall_distance(int n_, int*s1, int*s2){
    int dist = 0 ;
    int a, b;
    for(int i = 0 ; i < n_ ; i ++)
      for (int j=i+1; j < n_ ; j++){
         a =  s1[i] - s1[j];
         b =  s2[i] - s2[j];
         if (a*b < 0 ) dist +=1;
       }
    return dist;
}


int main(int argc, char const *argv[]) {
  /* code */
  int n_ = 8;
  int     *sigma_0 = new int[n_];
    int     *sigma_0_estim = new int[n_];
  kendall_distance(n_, sigma_0, sigma_0_estim);
  return 0;
}
