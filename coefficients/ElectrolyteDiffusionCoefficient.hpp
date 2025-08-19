#include "mfem.hpp"

using namespace mfem;

class DECoefficient : public Coefficient
{
private:
    const GridFunction &_u_gf;  
    double _ce_scale;             

public:
    DECoefficient(const GridFunction &u_gf, double ce_scale)
        : _u_gf(u_gf), _ce_scale(ce_scale) {}

    virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip) override
    {
<<<<<<< HEAD
        double u_val = _u_gf.GetValue(T, ip);
        //real_t u_val = 0.5;
        std::cout << "DE(u_val * _ce_scale) = " << DE(u_val * _ce_scale) << std::endl;

        return -DE(u_val * _ce_scale);
=======
        //double u_val = _u_gf.GetValue(T, ip);
        real_t u_val = 0.5;
        //std::cout << "u_val = " << u_val << std::endl;

        return DE(u_val * _ce_scale);
>>>>>>> 9cdeb5b93fb6d9d4942555f0815d3beba7cf25f8
    }
};
