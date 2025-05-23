namespace constants {
    enum Block : int {
        EP = 0,
        EC = 1,
        SP = 2,
        SC = 3
    };

    enum Region : int {
        PE  = 0,
        SEP = 1,
        NE  = 2
    };

    constexpr unsigned NPE = 5;
    constexpr unsigned NSEP = 5;
    constexpr unsigned NNE = 5;
    constexpr unsigned NX = NPE + NSEP + NNE;
    constexpr unsigned NR = 10;
    constexpr unsigned NPAR = (NPE - 1) + (NPE - 1);

    constexpr real_t LPE  = 1. * NPE  / NX;
    constexpr real_t LSEP = 1. * NSEP / NX;
    constexpr real_t LNE  = 1. * NNE  / NX;
}
