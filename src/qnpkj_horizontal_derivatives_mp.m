% DESCRIPTION: Computes horizontal derivatives of truncation coefficients
%              Q_{np}^{j} (Bucha et al., 2019a).  The Q_{np}^{j} coefficients
%              themselves can be evaluated using the
%              "qnpkj_radial_derivatives_mp" routine.
%
%              The code relies on the closed-form relations from Appendices
%              A.2, A.4, B.2 and B.4 of Bucha et al. (2019b) and uses the
%              outputs from "qnpkj_radial_derivatives_mp".
%
%              For an accurate evaluation, it is highly recommended to use an
%              extended number of significant digits when compared with double
%              precision. For instance, Bucha et al. (2019b) used 256
%              significant digits for n=0,...,21600 (harmonic degrees);
%              p=1,...,30 (topography integer power); and k=0,...,40 (order of
%              the radial derivative).
%
%              In Matlab, the number of significant digits can be extended via
%              the ADVANPIX Multiprecision Computing Toolbox
%              (https://www.advanpix.com/).  Full functional 7-day free trial
%              version of the toolbox for Microsoft Windows can be downloaded
%              from the web site.
%
%              This function requires the ADVANPIX toolbox. Nevertheless, after
%              some simple modifications, it could be rewritten into the
%              standard Matlab's syntax. This, however, may lead to inaccurate
%              results, depending on the variables "nmax", "pmax" and "kmax"
%              (see below, cf. also the Bucha et al. 2019a,b papers).
%
%
% INPUTS: Input parameters can be specified in the "Inputs" section below.
%        "nmax"    -- The maximum harmonic degree of the output truncation
%                     coefficients (nmax>=0)
%
%        "R"       -- Radius of the reference sphere to which the topography
%                     refers
%
%        "r"       -- Spherical radius of the evaluation point
%
%        "pmax"    -- Maximum number of the topography power (pmax>=1)
%
%        "psi"     -- Spherical distance separating the near- and far-zone 
%                     masses (inside and outside the spherical cap, respectively)
%
%        "kmax"    -- Maximum order of the radial derivative of the
%                     truncation coefficients (kmax>=0)
%
%        "near_zone_trunc_coeffs_path", "far_zone_trunc_coeffs_path" --
%                     Absolute or relative path to the near- and far-zone
%                     truncation coefficients from the
%                     "qnpkj_radial_derivatives_mp" function.  If you wish to
%                     compute horizontal derivatives for one type of
%                     coefficients only (either near- or far-zone ones), leave
%                     the other variable empty (e.g.,
%                     "far_zone_trunc_coeffs_path='';" will compute only
%                     near-zone coefficients)
%
%        "ADVANPIX_toolbox_path" -- Absolute or relative path to the ADVANPIX
%                     toolbox
%
%        "digits"  -- The number of significant digits used in the evaluation
%
%
% OUTPUTS: "save_name_path" -- Absolute or relative path to save the output
%                     near- and/or far-zone truncation coefficients files. The
%                     suffix ".mat" must not be included, because it is added
%                     automatically.

%                     The code generates and saves truncation coefficients
%                     "Q10", "Q11", "Q20", "Q21", "Q22", where "Qhi" stands for
%                     the h-th order horizontal derivative of the potential and
%                     "i" denotes derivatives with respect to the spherical
%                     distance "psi". For instance, "Q21" stands for truncation
%                     coefficients related to the gravity tensor (the
%                     second-order derivative of the potential) and the
%                     first-order derivative with respect to "psi".  The code
%                     automatically adds suffixes to the "save_name_path",
%                     depending on the particular truncation coefficients and
%                     whether these are related to near- or far-zone
%                     coefficients. For instance, far-zone coefficients "Q21"
%                     are saved in a binary mat-file with the name
%                     "save_name_path_Q21_far.mat"
%
%                     Depending on the inputs, the output variables may be
%                     vectors, 2D matrices or 3D matrices of dimensions
%                     (nmax+1,kmax+1,pmax). The first dimension provides
%                     truncation coefficients for harmonic degrees from 0 up to
%                     the user-defined "nmax" value. The second dimension
%                     yields radial derivatives of the truncation coefficients
%                     starting from the zero-order derivative up to the "kmax"
%                     order. The third dimension gives the coefficients for
%                     topography integer powers from 1 up to the user-defined
%                     "pmax" value.
%
%                     For instance, using "Q21" as an examples together with
%                     "nmax=2", "kmax=4" and "pmax=3" yields a 3D matrix with 
%                     the structure
%
%                     Q21(:,:,1) = [(0,0) (0,1) (0,2) (0,3) (0,4)]
%                                  [(1,0) (1,1) (1,2) (1,3) (1,4)]
%                                  [(2,0) (2,1) (2,2) (2,3) (2,4)]
%
%                     Q21(:,:,2) = [(0,0) (0,1) (0,2) (0,3) (0,4)]
%                                  [(1,0) (1,1) (1,2) (1,3) (1,4)]
%                                  [(2,0) (2,1) (2,2) (2,3) (2,4)]
%
%                     Q21(:,:,3) = [(0,0) (0,1) (0,2) (0,3) (0,4)]
%                                  [(1,0) (1,1) (1,2) (1,3) (1,4)]
%                                  [(2,0) (2,1) (2,2) (2,3) (2,4)]
%
%                     where (n,k) stands for harmonic degree "n" and order of
%                     the radial derivative "k". Note that a single layer of
%                     the "Q21" matrix (the third dimension) contains
%                     coefficients related to a fixed value of the topography
%                     integer power "p" (e.g., the third layer in the
%                     aforementioned example is related to coefficients for the
%                     third topography power "p").
%
% TEST RUN: Say we want to compute near- and far-zone truncation coefficients
%           "Q10", "Q11", "Q20", "Q21" and "Q22" up to degree 10, topography
%           power 3, order of the radial derivative 5, for a reference sphere
%           with the radius 6378137 m and at an evaluation point with the
%           radius 6378137 m + 7000 m. The integration radius of the near-zone
%           is defined by 100 km from the evaluation point.  During the
%           computation, we want to use 64 significant digits.  At first, we
%           need to evaluate and save the truncation coefficients "Q" using the
%           "qnpkj_radial_derivatives_mp" function. To this end, type following
%           commands in the Command Window:
%
%          addpath('/usr/local/MATLAB/Multiprecision Computing Toolbox/') %Path to the ADVANPIX toolbox
%          ndigits=64; %Defines the number of digits used in the computation
%          mp.Digits(ndigits); %Sets the number of significant digits to "ndigits"
%          psi=mp('100000')/mp('6378137'); %Radius of the integration zone
%          Q=qnpkj_radial_derivatives_mp(10,mp('6378137'),mp('6378137')+mp('7000'),3,cos(psi),7,0,ndigits);
%          mkdir('../data/test-outputs') %Creates output directory
%          save ../data/test-outputs/Truncation_coefficients_Q_near.mat -mat -v7.3 Q %saves the output Q coefficients
%          Q=qnpkj_radial_derivatives_mp(10,mp('6378137'),mp('6378137')+mp('7000'),3,cos(psi),7,1,ndigits);
%          save ../data/test-outputs/Truncation_coefficients_Q_far.mat -mat -v7.3 Q %saves the output Q coefficients
%
%          Note that to compute, for instance, the kth radial derivative of
%          "Q20" (in our test computation the 5th derivative), the coefficients 
%          from the "qnpkj_radial_derivatives_mp" function have to be computed 
%          up to the derivative "k+2" (here 7). 
%
%          Then, run this script while using the following input
%          parameters:

%               nmax   = 10;
%               R      = mp('6378137');
%               r      = R + mp('7000');
%               pmax   = 3;
%               psi    = mp('100000') / mp('6378137');
%               kmax   = 5;
%               digits = 64;

%          These values can be specified in the "Inputs" section below. Please
%          do not forget to set the correct path to the ADVANPIX toolbox.
%
%          The obtained values can be compared with the attached sample outputs
%          from the folder "../data/sample-outputs".
%
% REFERENCES: Bucha, B., Hirt, C., Kuhn, M. (2019a) Cap integration in spectral
%               gravity forward modelling: near- and far-zone gravity effects
%               via Molodensky's truncation coefficients. Journal of
%               Geodesy 93:65--83.
%
%            Bucha, B., Hirt, C., Kuhn, M. (2019b) Cap integration in
%               spectral gravity forward modelling up to the full gravity
%               tensor. Journal of Geodesy 93:1707--1737.
%
% CONTACT: blazej.bucha@stuba.sk
%
%
% Please use the following reference when using this function:
%
%         Bucha, B., Hirt, C., Kuhn, M., 2019. Cap integration in
%            spectral gravity forward modelling up to the full gravity
%            tensor. Journal of Geodesy 93:1707--1737.
%==========================================================================


clear; clc;


%Inputs
%==========================================================================
%Absolute or relative path to the ADVANPIX toolbox
ADVANPIX_toolbox_path='/usr/local/MATLAB/Multiprecision Computing Toolbox/';

digits=64; %The number of significant digits used in the computations. It
            %has to be the same as the total number of significant digits, with
            %which the beforehand prepared truncation coefficients "Q" and
            %their radial derivatives are provided

%No need to modify the next two commands
%..........................................................................
mp.Digits(digits); %The number of significant digits have to be set prior to any computations
addpath(ADVANPIX_toolbox_path);
%..........................................................................

nmax=10; %Maximum harmonic degree of truncation coefficients 
         %(must be equal to or smaller than that of the loaded
         %coefficients).
         %This value is an integer, so can be defined in the usual way,
         %that is, without the "mp" command.

R=mp('6378137'); %Radius of the reference sphere to which the topography refers
                 %(must be equal to that of the loaded coefficients).
                 %This value is generally not an integer, so have to be
                 %defined in the extended precision.

r=R+mp('7000'); %Spherical radius of the evaluation point
                %(must be equal to that of the loaded coefficients).
                %This value is generally not an integer, so have to be
                %defined in the extended precision.

pmax=3; %Maximum number of the topography power of the
        %output truncation coefficients (must be equal to or smaller 
        %than that of the loaded coefficients).
        %This value is an integer, so can be defined in the usual way,
        %that is, without the "mp" command.

psi=mp('100000')/mp('6378137'); %Spherical distance separating the near- and far-zone 
                                %masses (must be equal to that of the loaded coefficients)
                                %This value is generally not an integer, so have to be
                                %defined in the extended precision.

kmax=5; %Maximum order of the radial derivative of the truncation coefficients
        %(must be equal to that of the loaded coefficients).
        %This value is an integer, so can be defined in the usual way,
        %that is, without the "mp" command.

%Absolute or relative path to the near-zone truncation coefficients from
%the "qnpkj_radial_derivatives_mp" function. If you wish to compute
%horizontal derivatives only for the far-zone coefficients, leave this
%variable empty ("near_zone_trunc_coeffs_path='';")
near_zone_trunc_coeffs_path='../data/test-outputs/Truncation_coefficients_Q_near.mat';

%Absolute or relative path to the far-zone truncation coefficients from
%the "qnpkj_radial_derivatives_mp" function. If you wish to compute
%horizontal derivatives only for the near-zone coefficients, leave this
%variable empty ("far_zone_trunc_coeffs_path='';")
far_zone_trunc_coeffs_path='../data/test-outputs/Truncation_coefficients_Q_far.mat';
%==========================================================================


%Outputs
%==========================================================================
save_name_path='../data/test-outputs/Truncation_coefficients';
%==========================================================================


%Load truncation coefficients obtained by the "qnpkj_radial_derivatives_mp"
%function
%==========================================================================
if ~isempty(near_zone_trunc_coeffs_path)
    fprintf('Loading near-zone truncation coefficients... (%s) \n',datestr(clock));

    near_zone=1;
    c_near=1;

    Qn_near=load(near_zone_trunc_coeffs_path);
    Qn_near=struct2cell(Qn_near);
    Qn_near=cell2mat(Qn_near);
else
    near_zone=0;
end
if ~isempty(far_zone_trunc_coeffs_path)
    fprintf('Loading far-zone truncation coefficients... (%s) \n',datestr(clock));

    far_zone=1;
    c_far=-1;

    Qn_far=load(far_zone_trunc_coeffs_path);
    Qn_far=struct2cell(Qn_far);
    Qn_far=cell2mat(Qn_far);
else
    far_zone=0;
end
%==========================================================================


%Un-normalized associated Legendre functions
%==========================================================================
fprintf('Computing Legendre functions... (%s) \n',datestr(clock));

LF_order_max=2; %Maximum order of Legendre functions that will be evaluated
                      %The value is set to 2, because the code computes
                      %second-order horizontal derivatives

%Initialization
%--------------------------------------------------------------------------
cospsi=cos(psi);
sinpsi=sin(psi);
sqrt_1_m_cospsi2=sqrt(1-cospsi^2);
Pnm=zeros(nmax+1,LF_order_max+1,'mp');

P0=1;
P1=cospsi;
Pnm(1,1)=P0;
Pnm(2,1)=P1;
Pnm(2,2)=sinpsi;
%--------------------------------------------------------------------------


%Zonal Legendre functions
%--------------------------------------------------------------------------
for n=2:nmax
    P2=(2*n-1)/n*cospsi*P1-(n-1)/n*P0;
    Pnm(n+1,1)=P2;
    P0=P1;
    P1=P2;
end
%--------------------------------------------------------------------------


%Sectorial and tesseral Legendre functions
%--------------------------------------------------------------------------
for n=1:nmax
    for m=1:LF_order_max
        Pnm(n+1,m+1)=(Pnm(n,m+1)+(n-m+1)*sqrt_1_m_cospsi2*Pnm(n+1,m))/cospsi;
    end
end
%--------------------------------------------------------------------------
%==========================================================================


%Intergal kernels "Kp" and their radial derivatives via closed-form relations
%==========================================================================
max_der=pmax+kmax+2;


%Factorials and double factorials
%--------------------------------------------------------------------------
fprintf('Computing factorials and double factorial... (%s) \n',datestr(clock));
fact=factorial(0:(max_der+1));
double_fact=zeros(1,2*max_der,'mp');
for i=0:(2*max_der)
    double_fact(i+1)=prod(mp(i:-2:1));
end
%--------------------------------------------------------------------------


%Distance "l" and its radial derivatives (Eq. 68 of Bucha et al., 2019b)
%--------------------------------------------------------------------------
fprintf('Computing Euclidean distance and its derivatives... (%s) \n',datestr(clock));
l=zeros(1,max_der+1,'mp');
dist=sqrt(r^2-2*R*r*cospsi+R^2);
l(1)=1./dist;
for k=1:max_der
    for t=0:k
        if rem((k+t),2)==0
            l(k+1)=l(k+1)+(-1)^((k+t)/2)*double_fact(k-t+2)*double_fact(k+t)/fact(k-t+2)*fact(k+1)/fact(t+1)*((r-R*cospsi)^t)/dist^(k+t+1);
        end
    end
end
%--------------------------------------------------------------------------


%R-terms (Eq. 69 of Bucha et al., 2019b)
%--------------------------------------------------------------------------
fprintf('Computing the Rwq matrix... (%s) \n',datestr(clock));
Rwq=zeros(pmax,kmax+2,'mp');
for w=1:pmax
    for q=0:(kmax+1)
        if q==0
            rwq=r^w;
        else
            rwq=r.^(w-q);
            for j=1:q
                rwq=(w-j+1)*rwq;
            end
        end
        Rwq(w,q+1)=rwq;
    end
end
%--------------------------------------------------------------------------


%aps coefficients (Eq. 32 of Bucha et al., 2019b)
%--------------------------------------------------------------------------
fprintf('Computing the aps coefficients... (%s) \n',datestr(clock));
aps=zeros(pmax,pmax-2,'mp');
for p=1:pmax
    for s=1:p-2
        aps(p,s)=(-1)^(p-1)*fact(p)*fact(p-2)/fact(p-s+1)*fact(p-s-1)*fact(s);
    end
end
%--------------------------------------------------------------------------


%Binomial coefficients via recurrence relations using extended number of
%significant digits
%--------------------------------------------------------------------------
fprintf('Computing binomial coefficients... (%s) \n',datestr(clock));
binomial=zeros(max_der+1,max_der+1,'mp');
binomial(:,1)=1;
for i=0:max_der
    for ii=1:i
        binomial(i+1,ii+1)=binomial(i+1,ii)*(i-ii+1)/(ii);
    end
end
%--------------------------------------------------------------------------


%Integral kernels "Kp" and their radial derivatives (Eq. 67 of Bucha et al., 2019b)
%--------------------------------------------------------------------------
fprintf('Computing integration kernels using the closed relations... (%s) \n',datestr(clock));
K=zeros(pmax,kmax+2,'mp');
K1=R.*l;
K(1,1:(kmax+2))=K1(1,1:(kmax+2)); %Integral kernel for p=1
for k=0:(kmax+1) %Integral kernel for p=2
    K(2,k+1)=1/2*(-(k-1)*K1(1,k+1)-r*K1(1,k+2));
end
for k=0:(kmax+1)
    q=0:k;
    for p=3:pmax %Integral kernels for p>=3
        temp=0;
        for s=1:(p-2)
            temp=temp+aps(p,s)*sum(binomial(k+1,q+1).*Rwq(p-s,k-q+1).*K1(p-s+q+1));
        end
        K(p,k+1)=1/fact(p+1)*temp;
    end
end
%--------------------------------------------------------------------------
%==========================================================================


%Truncation coefficients for the first-order horizontal derivatives of the
%gravitational potential (gravitational vector)
%==========================================================================
%Truncation coefficients Q10
%--------------------------------------------------------------------------
fprintf('Computing truncation coefficients Q10... (%s) \n',datestr(clock));

%(Eq. 65 of Bucha et al., 2019b)
if near_zone==1
    Q10_near=Qn_near(1:(nmax+1),2:(kmax+1),:); %#ok<NASGU>
    eval(sprintf('save %s_Q10_near.mat -mat -v7.3 Q10_near',save_name_path))
    clear Q10_near
end
if far_zone==1
    Q10_far=Qn_far(1:(nmax+1),2:(kmax+1),:); %#ok<NASGU>
    eval(sprintf('save %s_Q10_far.mat -mat -v7.3 Q10_far',save_name_path))
    clear Q10_far
end
%--------------------------------------------------------------------------


%Truncation coeffcients Q11
%--------------------------------------------------------------------------
fprintf('Computing truncation coefficients Q11... (%s) \n',datestr(clock));

if near_zone==1
    Q11_near=zeros(nmax+1,kmax+1,pmax,'mp');
end
if far_zone==1
    Q11_far=zeros(nmax+1,kmax+1,pmax,'mp');
end
n=0:nmax;
n=n(:);
sinpsi_Pn1_nnp1=sinpsi.*Pnm(:,2)./(n.*(n+1));

%(Eq. 66 of Bucha et al., 2019b)
for p=1:pmax
    for k=0:kmax
        temp_near=0;
        temp_far=0;
        for q=0:k
            temp=binomial(k+1,q+1)*(-1)^(k-q)*fact(k-q+1)./r^(k-q+1);
            temp2=sinpsi_Pn1_nnp1.*K(p,q+1);
            if near_zone==1
                temp_near=temp_near+temp.*(c_near*temp2-Qn_near(1:(nmax+1),q+1,p));
            end
            if far_zone==1
                temp_far=temp_far+temp.*(c_far*temp2-Qn_far(1:(nmax+1),q+1,p));
            end
        end
        Q11_near(1:(nmax+1),k+1,p)=temp_near;
        Q11_far(1:(nmax+1),k+1,p)=temp_far;
    end
end
clear temp temp2 temp_near temp_far
if near_zone==1
    Q11_near(1,:,:)=0; %These coefficients do not exist, so are here set to zero
    eval(sprintf('save %s_Q11_near.mat -mat -v7.3 Q11_near',save_name_path))
end
if far_zone==1
    Q11_far(1,:,:)=0; %These coefficients do not exist, so are here set to zero
    eval(sprintf('save %s_Q11_far.mat -mat -v7.3 Q11_far',save_name_path))
end
%--------------------------------------------------------------------------
%==========================================================================


%Truncation coefficients for the secondorder horizontal derivatives of the
%gravitational potential (gravitational tensor)
%==========================================================================
%Truncation coefficients Q20
%--------------------------------------------------------------------------
fprintf('Computing truncation coefficients Q20... (%s) \n',datestr(clock));

%(Eq. 97 of Bucha et al., 2019b)
if near_zone==1
    Q20_near=Qn_near(1:(nmax+1),3:(kmax+3),:); %#ok<NASGU>
    eval(sprintf('save %s_Q20_near.mat -mat -v7.3 Q20_near',save_name_path))
    clear Q20_near
end
if far_zone==1
    Q20_far=Qn_far(1:(nmax+1),3:(kmax+3),:); %#ok<NASGU>
    eval(sprintf('save %s_Q20_far.mat -mat -v7.3 Q20_far',save_name_path))
    clear Q20_far
end
%--------------------------------------------------------------------------


%Truncation coefficients Q21
%--------------------------------------------------------------------------
fprintf('Computing truncation coefficients Q21... (%s) \n',datestr(clock));

if near_zone==1
    Q21_near=zeros(nmax+1,kmax+1,pmax,'mp');
end
if far_zone==1
    Q21_far=zeros(nmax+1,kmax+1,pmax,'mp');
end

%(Eq. 98 of Bucha et al., 2019b)
for p=1:pmax
    for k=0:kmax
        temp_near=0;
        temp_far=0;
        for q=0:k
            temp=binomial(k+1,q+1)*(-1)^(k-q)*fact(k-q+1)./r^(k-q+1);
            temp2=sinpsi_Pn1_nnp1.*K(p,q+2);
            if near_zone==1
                temp_near=temp_near+temp.*(Q11_near(1:(nmax+1),q+1,p)-c_near*temp2+Qn_near(1:(nmax+1),q+2,p));
            end
            if far_zone==1
                temp_far=temp_far+temp.*(Q11_far(1:(nmax+1),q+1,p)-c_far*temp2+Qn_far(1:(nmax+1),q+2,p));
            end
        end
        Q21_near(1:(nmax+1),k+1,p)=temp_near;
        Q21_far(1:(nmax+1),k+1,p)=temp_far;
    end
end
clear temp temp2 temp_near temp_far
if near_zone==1
    clear Q11_near
    Q21_near(1,:,:)=0; %#ok<NASGU> %These coefficients do not exist, so are here set to zero
    eval(sprintf('save %s_Q21_near.mat -mat -v7.3 Q21_near',save_name_path))
    clear Q21_near
end
if far_zone==1
    clear Q11_far
    Q21_far(1,:,:)=0; %#ok<NASGU> %These coefficients do not exist, so are here set to zero
    eval(sprintf('save %s_Q21_far.mat -mat -v7.3 Q21_far',save_name_path))
    clear Q21_far
end
%--------------------------------------------------------------------------


%Truncation coefficients Q22
%--------------------------------------------------------------------------
fprintf('Computing truncation coefficients Q22... (%s) \n',datestr(clock));

%Derivatives of the distance "l" with respect to "r" (multiple times) and
%"psi" (only once) (Eq. 101 of Bucha et al., 2019b)
%..........................................................................
lpsi=zeros(1,max_der+1,'mp');
lpsi(1)=-(R*r*sinpsi)/dist^3;
for k=1:max_der
    for t=0:k
        if rem((k+t),2)==0
            lpsi(k+1)=lpsi(k+1)+(-1)^((k+t)/2)*double_fact(k-t+2)*double_fact(k+t)/fact(k-t+2)*fact(k+1)/fact(t+1)*(((r-R*cospsi)^(t-1)*R*sinpsi)/dist^(k+t+1)*(t-(r-R*cospsi)*r*(k+t+1)/(dist^2)));
        end
    end
end
%..........................................................................


%Derivatives of the integral kernels "Kp" with respect to "r" (multiple times) and
%"psi" (only once) (Eq. 100 of Bucha et al., 2019b)
%..........................................................................
Kpsi=zeros(pmax,kmax+2,'mp');
Kpsi1=R*lpsi;
Kpsi(1,1:(kmax+2))=Kpsi1(1,1:(kmax+2)); %p=1
for k=0:(kmax+1) %p=2
    Kpsi(2,k+1)=1/2*(-(k-1)*Kpsi1(1,k+1)-r*Kpsi1(1,k+2));
end
for k=0:(kmax+1)
    q=0:k;
    for p=3:pmax %p>=3
        temp=0;
        for s=1:(p-2)
            temp=temp+aps(p,s)*sum(binomial(k+1,q+1).*Rwq(p-s,k-q+1).*Kpsi1(p-s+q+1));
        end
        Kpsi(p,k+1)=1/fact(p+1)*temp;
    end
end
%..........................................................................


if near_zone==1
    Q22_near=zeros(nmax+1,kmax+1,pmax,'mp');
end
if far_zone==1
    Q22_far=zeros(nmax+1,kmax+1,pmax,'mp');
end
sinpsi_Pn2_np2p1nm1=sinpsi.*Pnm(:,3)./((n+2).*(n+1).*n.*(n-1));

%(Eq. 99 of Bucha et al., 2019b)
for p=1:pmax
    for k=0:kmax
        temp_near=0;
        temp_far=0;
        for q=0:k
            temp=binomial(k+1,q+1)*(-1)^(k-q)*fact(k-q+2)./r^(k-q+2);
            temp2=sinpsi_Pn2_np2p1nm1.*Kpsi(p,q+1);
            temp3=sinpsi_Pn1_nnp1.*K(p,q+1);
            if near_zone==1
                temp_near=temp_near+temp.*(c_near*temp2-c_near*temp3+Qn_near(1:(nmax+1),q+1,p));
            end
            if far_zone==1
                temp_far=temp_far+temp.*(c_far*temp2-c_far*temp3+Qn_far(1:(nmax+1),q+1,p));
            end
        end
        Q22_near(1:(nmax+1),k+1,p)=1/2*temp_near;
        Q22_far(1:(nmax+1),k+1,p)=1/2*temp_far;
    end
end
clear temp temp2 temp3 temp_near temp_far
if near_zone==1
    Q22_near(1:2,:,:)=0; %#ok<NASGU> %These coefficients do not exist, so are here set to zero
    eval(sprintf('save %s_Q22_near.mat -mat -v7.3 Q22_near',save_name_path))
    clear Q22_near
end
if far_zone==1
    Q22_far(1:2,:,:)=0; %#ok<NASGU> %These coefficients do not exist, so are here set to zero
    eval(sprintf('save %s_Q22_far.mat -mat -v7.3 Q22_far',save_name_path))
    clear Q22_far
end
%--------------------------------------------------------------------------
%==========================================================================

fprintf('The computation has been finished... (%s) \n',datestr(clock));
