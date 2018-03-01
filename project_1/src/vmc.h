class vmc {
    int n_dim, N, n_carlos;
    double step;

    public :
        vmc(int n_dim, int N, int n_carlos, double step);
        void solve(string sampling, string derivative, bool variatonal_search, bool interact);
        void output_result();

    private :
        carlo_cycle(mat R, mat R_p, double alpha, double beta, bool importance);
        variational_mc(double beta_start, double beta_stop, double beta_increment,
                double alpha_start, double alpha_stop, double alpha_increment, bool gradient_desc);

};
