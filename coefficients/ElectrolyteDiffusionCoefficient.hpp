#include "mfem.hpp"

using namespace mfem;

class DECoefficient : public Coefficient
{
private:
    const GridFunction &_u_gf;  
    double _ce_scale;             

public:
    DECoefficient(const GridFunction &u_gf)
        : _u_gf(u_gf) {}

    virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip) override
    {
        double u_val = _u_gf.GetValue(T, ip) + CE0;
        return DE(u_val);
    }
};
