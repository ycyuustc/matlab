%% Grassmann Numbers and Clifford Algebras
% _The elementary particles of Physics are classified according to the behavior
% of the multi-particle states under exchange of identical particles: bosonic
% states are symmetric while fermionic This manifests itself also in the commutation
% properties of the respective creation operators: bosonic creation operators
% commute while fermionic ones anticommute._
%
% The function fn_c represents the function fn_contact which is introduced
% in the MPS calculations, and the function fn_f are the function in spin theories
% which is used to generate the $<math xmlns="http://www.w3.org/1998/Math/MathML"
% display="inline"><mrow><msup><mrow><mi mathvariant="italic">s</mi></mrow><mrow><mo>+</mo></mrow></msup></mrow></math>$
% and $<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mrow><msup><mrow><mi
% mathvariant="italic">s</mi></mrow><mrow><mo>?</mo></mrow></msup></mrow></math>$
% presentations.

fn_f=@(j,m) sqrt((j+m).*(j-m+1));
fn_c=@(x1,x2,x3,x4,x5,x6) fn_contract(x1,x2,x3,x4,x5,x6);
%%
% the total spin J=1/2; and Num_J is the dimension of the representation
% for each site.
%%
J=1/2;    Num_J=2*J+1;
%%
%  then we generate the Z operator:

v_J=(Num_J-1)/2:-1:-(Num_J-1)/2; mJ_z=diag(v_J);
%%
% and then we can introduce the $<math xmlns="http://www.w3.org/1998/Math/MathML"
% display="inline"><mrow><msup><mrow><mi mathvariant="italic">J</mi></mrow><mrow><mo>+</mo></mrow></msup><mo>,</mo><msup><mrow><mi
% mathvariant="italic">J</mi></mrow><mrow><mo>?</mo></mrow></msup><mo>,</mo><msub><mrow><mi
% mathvariant="italic">J</mi></mrow><mrow><mi mathvariant="italic">x</mi></mrow></msub></mrow></math>$
% and $<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mrow><msub><mrow><mi
% mathvariant="italic">J</mi></mrow><mrow><mi mathvariant="italic">y</mi></mrow></msub></mrow></math>$
% operator

v_m=J:-1:-J+1;
mJ_p=diag(fn_f(J,v_m),1); mJ_m=diag(fn_f(J,v_m),-1);
mJ_x=(mJ_p+mJ_m)/2; mJ_y=(mJ_p-mJ_m)/(2j);
%%
% Finally we get the identity operator $<math xmlns="http://www.w3.org/1998/Math/MathML"
% display="inline"><mrow><mi mathvariant="italic">I</mi></mrow></math>$ :

mJ_id=eye(Num_J);
%% The delta matrix
% The delta matrices are closely related to the spin operator, we wish we could
% have:
%
% $$[ \sigma_i, \sigma_j ] = 2 \epsilon_{ijk} \sigma_k$$
%%
m_dtx=mJ_x*2;
m_dty=mJ_y*2;
m_dtz=mJ_z*2;
m_dti=mJ_id;
%% Show the above result:
% mJ_p is:
%%
% disp('mJ_p=');
disp(mJ_p); disp(2); disp(3)
%%
% mJ_m is:

disp('mJ_m=');disp(mJ_m);
%%
% identity operator mJ_id is:

disp('mJ_id');disp(mJ_id);
%%
% The delta matrices are:

m_dtx
m_dty
m_dtz
m_dti
%%