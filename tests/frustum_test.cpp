#include "frustum_test.hpp"

#include <vmmlib/frustum.hpp>
#include <vmmlib/frustum_culler.hpp>

namespace vmml
{

bool
frustum_test::run() 
{
    {
        // TODO
        frustum< float >        f;
        frustum_culler< float > fc;
        
        matrix< 4, 4, float >   mvp;
        
        fc.setup( mvp );
    
    }

    return true;
}


} // namespace vmml

