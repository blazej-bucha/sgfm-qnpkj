function Q=qnpkj_radial_derivatives_mp(nmax,R,r,pmax,cos_psi,kmax,zone,digits)
% ==========================================================================
% DESCRIPTION: Computes truncation coefficients Q_{np}^{j} and their radial
%              derivatives from the Bucha et al. (2019) paper via recurrence
%              relations with a fixed number of terms (Appendices C and D).
%
%              For an accurate evaluation, it is highly recommended to use an
%              extended number of significant digits when compared with double
%              precision.  For instance, Bucha et al. (2019) used 64
%              significant digits for n=0,...,21600 (harmonic degrees);
%              p=1,...,15 (topography integer power); and k=0,...,15 (order of
%              the radial derivative).  For higher values of "p" and "k", even
%              larger number of significant digits may need to be employed
%              (e.g., 256 when "pmax=30" and "kmax=40").
%
%              In Matlab, the number of significant digits can be extended via
%              the ADVANPIX Multiprecision Computing Toolbox
%              (https://www.advanpix.com/).  Full functional 7-day free trial
%              version of the toolbox for Microsoft Windows can be downloaded
%              from the web site.
%
%              This function requires the ADVANPIX toolbox.  Nevertheless,
%              after some simple modifications, it could be rewritten into the
%              standard Matlab's syntax.  This, however, may lead to inaccurate
%              results, depending on the variables "nmax", "pmax" and "kmax"
%              (see below, cf. also the Bucha et al. 2019 paper).
%
% INPUTS: "nmax"    -- The maximum harmonic degree of the output truncation
%                      coefficients (nmax>=0)
%
%         "R"       -- Radius of the reference sphere to which the topography
%                      refers
%
%         "r"       -- Spherical radius of the evaluation point
%
%         "pmax"    -- Maximum number of the topography power (pmax>=1)
%
%         "cos_psi" -- Cosine of the spherical distance separating the near-
%                      and far-zone masses (inside and outside the spherical
%                      cap, respectively)
%
%         "kmax"    -- Maximum order of the radial derivative of the truncation
%                      coefficients (kmax>=0)
%
%         "zone"    -- 0 - Truncation coefficients for near-zone effects
%                      1 - Truncation coefficients for far-zone effects
%
%         "digits"  -- The number of significant digits used in the evaluation
%
%
%OUTPUTS: "Q" -- Truncation coefficients.  Depending on the inputs, "Q"
%                may be a vector, 2D matrix or a 3D matrix of dimensions
%                (nmax+1,kmax+1,pmax).  The first dimension provides
%                truncation coefficients for harmonic degrees from 0 up to
%                the user-defined "nmax" value.  The second dimension yields
%                radial derivatives of the truncation coefficients starting
%                from the zero-order derivative up to the "kmax" order.  The
%                third dimension gives the coefficients for topography
%                integer powers from 1 up to the user-defined "pmax" value.
%
%                For instance, using "nmax=2", "kmax=4" and
%                "pmax=3" yields a 3D matrix with the structure
%
%                   Q(:,:,1) = [(0,0) (0,1) (0,2) (0,3) (0,4)]
%                              [(1,0) (1,1) (1,2) (1,3) (1,4)]
%                              [(2,0) (2,1) (2,2) (2,3) (2,4)]
%
%                   Q(:,:,2) = [(0,0) (0,1) (0,2) (0,3) (0,4)]
%                              [(1,0) (1,1) (1,2) (1,3) (1,4)]
%                              [(2,0) (2,1) (2,2) (2,3) (2,4)]
%
%                   Q(:,:,3) = [(0,0) (0,1) (0,2) (0,3) (0,4)]
%                              [(1,0) (1,1) (1,2) (1,3) (1,4)]
%                              [(2,0) (2,1) (2,2) (2,3) (2,4)]
%
%                where (n,k) stands for harmonic degree "n" and order of
%                the radial derivative "k".  Note that a single layer of the
%                "Q" matrix (the third dimension) contains coefficients
%                related to a fixed value of the topography integer power
%                "p" (e.g., the third layer in the aforementioned example is
%                related to coefficients for the third topography power "p").
%
% TEST RUN: Say we want to compute near-zone truncation coefficients up to
%           degree 10, topography power 3, order of the radial derivative 7,
%           for a reference sphere with the radius 6378137 m and at an
%           evaluation point with the radius 6378137 m + 7000 m.  The
%           integration radius of the near-zone is defined by 100 km from the
%           evaluation point.  During the computation, we want to use 64
%           significant digits.  In the Command Window, type following five
%           commands:
%
%           addpath('/usr/local/MATLAB/Multiprecision Computing Toolbox/') %Path to the ADVANPIX toolbox
%           ndigits=64; %Defines the number of digits used in the computation
%           mp.Digits(ndigits); %Sets the number of significant digits to "ndigits"
%           psi=mp('100000')/mp('6378137'); %Radius of the integration zone
%           Q=qnpkj_radial_derivatives_mp(10,mp('6378137'),mp('6378137')+mp('7000'),3,cos(psi),7,0,ndigits);
%
%           The obtained values can be compared with the attached sample
%           outputs from the file
%           "sample-outputs/Truncation_coefficients_Q_near.mat".
%
%           Importantly, note that the number of significant digits needs to be
%           set also outside the "qnpkj_radial_derivatives_mp" function.  This
%           is because the input variable "cos_psi" is generally not an
%           integer, but instead a real number.  In our case, the spherical
%           distance "psi" reads in radians "psi=100000/6378137", which is not
%           an integer.  "psi" must therefore be prepared using an extended
%           number of significant digits as "mp('100000')/mp('6378137')", and
%           only after this it can be used as an input.  The same holds true
%           for the "R" and "r" variables.
%
%           Also note that integer input variables "nmax", "pmax", "kmax",
%           "zone" and "digits" may be defined in the usual way (double
%           precision), that is, without the "mp" command.  In fact, this is
%           the recommended way in order not to slow down the computation.
%
% REFERENCES: Bucha, B., Hirt, C., Kuhn, M. (2019) Cap integration in spectral
%               gravity forward modelling: near- and far-zone gravity effects
%               via Molodensky's truncation coefficients.  Journal of Geodesy
%               93:65--83.
%
% CONTACT: blazej.bucha@stuba.sk
%
%
% Please use the following reference when using this function:
%
%         Bucha, B., Hirt, C., Kuhn, M., 2019.  Cap integration in
%            spectral gravity forward modelling up to the full gravity
%            tensor.  Journal of Geodesy 93:1707--1737.
% ==========================================================================

mp.Digits(digits);

%Useful substitution
%--------------------------------------------------------------------------
mpnmax=mp(nmax);
mp1d3=mp('1/3');
%--------------------------------------------------------------------------


%Substitutions and initializations
%--------------------------------------------------------------------------
t=R/r; %Substitution
t2=t^2; %Substitution
max_der=pmax+kmax+1; %Order of the maximum derivative that needs to be 
                       %computed in some parts of the code

if zone==0 %Near zone
    ul=cos_psi; %Lower boundary of the integration over "psi"
    uu=1; %Upper boundary of the integration over "psi"
elseif zone==1 %Far zone
    ul=-1; %Lower boundary of the integration over "psi"
    uu=cos_psi; %Upper boundary of the integration over "psi"
end
gu_der=zeros(1,max_der+1,'mp'); %Initialization of a matrix with derivatives
                                %of "1/g(t,uu)" with respect to "t" (Eq. 61 of Bucha et al., 2019)
gl_der=zeros(1,max_der+1,'mp'); %Initialization of a matrix with derivatives
                                %of "1/g(t,ul)" with respect to "t" (Eq. 61 of Bucha et al., 2019)
%--------------------------------------------------------------------------


%Factorials
%--------------------------------------------------------------------------
fact=factorial(mp(0:(max_der+1)));
double_fact=zeros(1,2*max_der,'mp');
for i=0:(2*max_der)
    double_fact(i+1)=prod(mp(i:-2:1));
end
%--------------------------------------------------------------------------


%Derivatives of "1/g(t,uu)" with respect to "t" (Eq. 61 of Bucha et al., 2019)
%--------------------------------------------------------------------------
fprintf('Computing gu-terms (%s)\n',datestr(clock))

gu=sqrt(1-2*t*uu+t2); %Normalized distance "g(t,uu)"
gu_der(1)=1/gu; %Zero-order derivative of "1/g(t,uu)" with respect to "t"
for i=1:max_der %First- and higher-order derivative of "1/g(t,uu)" with respect to "t"
    for l=0:i
        if rem((i+l),2)==0
            gu_der(i+1)=gu_der(i+1)+(-1)^((i+l)/2)*double_fact(i-l+2)*double_fact(i+l)/fact(i-l+2)*fact(i+1)/fact(l+1)*((t-uu)^l)/gu^(i+l+1);
        end
    end
end
%--------------------------------------------------------------------------


%Derivatives of "1/g(t,ul)" with respect to "t" (Eq. 61 of Bucha et al., 2019)
%--------------------------------------------------------------------------
fprintf('Computing gl-terms (%s)\n',datestr(clock))

gl=sqrt(1-2*t*ul+t2); %Normalized distance "g(t,ul)"
gl_der(1)=1/gl; %Zero-order derivative of "1/g(t,ul)" with respect to "t"
for i=1:max_der %First- and higher-order derivative of "1/g(t,ul)" with respect to "t"
    for l=0:i
        if rem((i+l),2)==0
            gl_der(i+1)=gl_der(i+1)+(-1)^((i+l)/2)*double_fact(i-l+2)*double_fact(i+l)/fact(i-l+2)*fact(i+1)/fact(l+1)*((t-ul)^l)/gl^(i+l+1);
        end
    end
end
%--------------------------------------------------------------------------


if zone==0 %Near zone
    u0=ul;
    g0_der=gl_der;
elseif zone==1 %Far zone
    u0=uu;
    g0_der=gu_der;
end


%Integrals of Legendre polynomials (Eqs. 69 -- 71 of Bucha et al., 2019)
%--------------------------------------------------------------------------
fprintf('Computing M-terms (%s)\n',datestr(clock))

%Initializations
M=zeros(1,nmax+1,'mp');
if zone==0 %Near zone
    M(1)=1-ul;
    M(2)=0.5*(1-ul^2);
elseif zone==1 %Far zone
    M(1)=1+uu;
    M(2)=0.5*(uu^2-1);
end

nn=0:mpnmax;
nn1=1./(nn+1); %Substitution
nn2=2.*nn-1; %Substitution
nnm2=nn-2; %Substitution
for n=2:nmax %Loop over harmonic degrees
    M(n+1)=nn1(n+1)*(nn2(n+1)*u0*M(n)-nnm2(n+1)*M(n-1));
end
clear nn1 nn2 nnm2 %These terms are no longer necessary
%--------------------------------------------------------------------------


%Binomial coefficients
%--------------------------------------------------------------------------
fprintf('Computing binomial coefficients (%s)\n',datestr(clock))

%Binomial coefficients are here computed via recurrence relations using the
%extended number of significant digits
binomial=zeros(max_der+1,max_der+1,'mp');
binomial(:,1)=1;
for i=0:max_der
    for ii=1:i
        binomial(i+1,ii+1)=binomial(i+1,ii)*(i-ii+1)/ii;
    end
end
%--------------------------------------------------------------------------


%B-terms (Eqs. 64 -- 66 of Bucha et al., 2019)
%--------------------------------------------------------------------------
fprintf('Computing B-terms (%s)\n',datestr(clock))

%Initializations
B=zeros(nmax+1,max_der+1,'mp');
t_temp=(1+t2)/t; %Substitution
beta0uu=zeros(1,max_der+1,'mp');
beta1uu=beta0uu;
beta1ul=beta0uu;
for i=0:max_der
    tmi1=t^(-(i+1));
    tmi2=t^(-(i+2));
    m1i=(-1)^i;
    m1i1=(-1)^(i+1);
    tp2=1/t^2;
    beta0uu(i+1)=m1i*fact(i+1)*tmi1;
    if i==0
        beta1uu(i+1)=tp2-uu/t+1;
        beta1ul(i+1)=tp2-ul/t+1;
    else
        fact_i1=fact(i+1);
        fact_i2=fact(i+2);
        beta1uu(i+1)=m1i*fact_i2*tmi2+m1i1*fact_i1*tmi1*uu;
        beta1ul(i+1)=m1i*fact_i2*tmi2+m1i1*fact_i1*tmi1*ul;
    end
end
beta0ul=beta0uu;

iii=0:max_der; %Substitution
it=iii./t; %Substitution
iim1t=iii.*(iii-1)./t; %Substitution
iii2=2*iii; %Substitution
g0_dert=g0_der./t; %Substitution
for i=0:max_der %Loop over derivatives
    %Initial values
    if i==0
        B(1,i+1)=1/(t*gu)-1/(t*gl);
        if nmax>0
            B(2,i+1)=1/(t2*gu)*(1-t*uu+t2)-1/(t2*gl)*(1-t*ul+t2);
        end
    else
        for k=0:i
            ii=i-k;

            B(1,i+1)=B(1,i+1)+binomial(i+1,k+1)*(beta0uu(ii+1)*gu_der(k+1)-beta0ul(ii+1)*gl_der(k+1));
            if nmax>0
                B(2,i+1)=B(2,i+1)+binomial(i+1,k+1)*(beta1uu(ii+1)*gu_der(k+1)-beta1ul(ii+1)*gl_der(k+1));
            end
        end
    end

    it_i=it(i+1); %Substitution
    iii2_i=iii2(i+1); %Substitution
    iim1t_i=iim1t(i+1); %Substitution
    Mn_g0_dert_i=M*g0_dert(i+1); %Substitution
    for n=2:nmax %Loop over harmonic degrees
        if i==0
            B(n+1,i+1)=t_temp*B(n,i+1)-B(n-1,i+1)-Mn_g0_dert_i(n);
        elseif i==1
            B(n+1,i+1)=t_temp*B(n,i+1)-B(n-1,i+1)-it_i*(B(n+1,i)+B(n-1,i))+iii2_i*B(n,i)-Mn_g0_dert_i(n);
        else
            B(n+1,i+1)=t_temp*B(n,i+1)-B(n-1,i+1)-it_i*(B(n+1,i)+B(n-1,i))+iii2_i*B(n,i)+iim1t_i*B(n,i-1)-Mn_g0_dert_i(n);
        end
    end
end
%--------------------------------------------------------------------------


%A-terms (Eqs. 56 -- 58 of Bucha et al., 2019)
%--------------------------------------------------------------------------
fprintf('Computing A-terms (%s)\n',datestr(clock))

%Initializations
A=zeros(nmax+1,max_der+1,'mp');
alfa0uu=zeros(1,max_der+1,'mp');
alfa0ul=alfa0uu;
alfa1uu=alfa0uu;
alfa1ul=alfa0uu;
for i=0:max_der
    m1i=(-1)^i;
    m1i1=(-1)^(i+1);
    fact_i1=fact(i+1);
    fact_i2=fact(i+2);
    ti1=t^(i+1);
    ti2=t^(i+2);
    m1i_fact_i1_ti1=m1i*fact_i1/ti1;
    m1i_fact_i2_ti2=m1i*fact_i2/ti2;
    if i==0
        alfa0uu(i+1)=(1/t-2*uu+t);
        alfa0ul(i+1)=(1/t-2*ul+t);
        alfa1uu(i+1)=t2-t*uu-uu/t+1/t2-2*uu^2+2;
        alfa1ul(i+1)=t2-t*ul-ul/t+1/t2-2*ul^2+2;
    elseif i==1
        alfa0uu(i+1)=m1i_fact_i1_ti1+1;
        alfa0ul(i+1)=m1i_fact_i1_ti1+1;
        alfa1uu(i+1)=m1i1*fact_i1*uu/ti1+m1i_fact_i2_ti2+2*t-uu;
        alfa1ul(i+1)=m1i1*fact_i1*ul/ti1+m1i_fact_i2_ti2+2*t-ul;
    elseif i==2
        alfa0uu(i+1)=m1i_fact_i1_ti1;
        alfa0ul(i+1)=m1i_fact_i1_ti1;
        alfa1uu(i+1)=m1i1*fact_i1*uu/ti1+m1i_fact_i2_ti2+2;
        alfa1ul(i+1)=m1i1*fact_i1*ul/ti1+m1i_fact_i2_ti2+2;
    else
        alfa0uu(i+1)=m1i_fact_i1_ti1;
        alfa0ul(i+1)=m1i_fact_i1_ti1;
        alfa1uu(i+1)=m1i1*fact_i1*uu/ti1+m1i_fact_i2_ti2;
        alfa1ul(i+1)=m1i1*fact_i1*ul/ti1+m1i_fact_i2_ti2;
    end
end

nn=0:mpnmax; %Substitution
nn21=1./(2*nn(3:end)'+1); %Substitution
nn2=2*(nn(3:end)'+1); %Substitution
t2p1=1+t2; %Substitution
tt2=2*t; %Substitution
ii_1=iii.*(iii-1); %Substitution
M=M(3:end)'; %Substitution
nn2M=nn2.*M; %Substitution
for i=0:max_der
    if i==0
        A(1,i+1)=-gu/t+gl/t;
        if nmax>0
            A(2,i+1)=-gu/(3*t2)*(1+t*uu+t2)+gl/(3*t2)*(1+t*ul+t2);
        end
    else
        for k=0:i
            ii=i-k;

            A(1,i+1)=A(1,i+1)+binomial(i+1,k+1)*(alfa0uu(ii+1)*gu_der(k+1)-alfa0ul(ii+1)*gl_der(k+1));
            if nmax>0
                A(2,i+1)=A(2,i+1)+binomial(i+1,k+1)*(alfa1uu(ii+1)*gu_der(k+1)-alfa1ul(ii+1)*gl_der(k+1));
            end
        end
        A(1,i+1)=-A(1,i+1);
        if nmax>0
            A(2,i+1)=-mp1d3*A(2,i+1);
        end
    end

    Bnp1_ip1=B(3:end,i+1);
    Bn_ip1=B(2:end-1,i+1);
    if i>=2
        Bnp1_im1=B(3:end,i-1);
    end
    if i>=1
        Bnp1_i=B(3:end,i);
        Bn_i=B(2:end-1,i);
    end
    if i==0
        A(3:end,i+1)=nn21.*(-(t2p1).*Bnp1_ip1+tt2.*Bn_ip1+nn2M*g0_der(i+1));
    elseif i==1
        A(3:end,i+1)=nn21.*(-(t2p1).*Bnp1_ip1+tt2.*(Bn_ip1-i*Bnp1_i)+iii2(i+1).*Bn_i+nn2M*g0_der(i+1));
    else
        A(3:end,i+1)=nn21.*(-(t2p1).*Bnp1_ip1+tt2.*(Bn_ip1-i*Bnp1_i)+iii2(i+1).*Bn_i-ii_1(i+1).*Bnp1_im1+nn2M*g0_der(i+1));
    end
end
clear B Bn_i Bn_ip1 Bnp1_im1 Bnp1_ip1 M Mn_g0_dert_i nn nn21 nn2 t2p1 tt2 iii2 ii_1 M nn2M %These are no longer necessary
%--------------------------------------------------------------------------


%L-terms (Eq. 46 of Bucha et al., 2019)
%--------------------------------------------------------------------------
fprintf('Computing L-terms (%s)\n',datestr(clock))

L=zeros(nmax+1,max_der+1,'mp');
for i=0:max_der
    if i==0
        L(:,i+1)=t*A(:,i+1);
    else
        L(:,i+1)=t*A(:,i+1)+i*A(:,i);
    end
end
%--------------------------------------------------------------------------


%T-terms (Eq. 50 of Bucha et al., 2019)
%--------------------------------------------------------------------------
fprintf('Computing T-terms (%s)\n',datestr(clock))

T=zeros(max_der+1,1,'mp');
for i=0:max_der
    T(i+1)=(-1)^i*fact(i+1)*R/r^(i+1);
end
%--------------------------------------------------------------------------


%G-terms (Eq. 48 of Bucha et al., 2019)
%--------------------------------------------------------------------------
fprintf('Computing G-terms (%s)\n',datestr(clock))

G=zeros(nmax+1,max_der+1,'mp');
for i=0:max_der
    if i==0
        G(:,i+1)=t*A(:,i+1);
    else
        %Partial Bell polynomials
        Bell=zeros(max_der+1,max_der+1,'mp');
        Bell(1,1)=1;
        Bell(2:end,1)=0;
        Bell(1,2:end)=0;
        for n=1:max_der
            for k=1:max_der
                if n-k+1>0
                    Bell(n+1,k+1)=Bell(n+1,k+1)+(binomial(n,1:(n-k+1))*(T((1:(n-k+1))+1).*Bell(n-(1:(n-k+1))+1,k)));
                end
            end
        end

        for pp=1:i
            G(:,i+1)=G(:,i+1)+L(:,pp+1)*Bell(i+1,pp+1);
        end
    end
end
clear A T Bell %These terms are no longer necessary
%--------------------------------------------------------------------------


%Truncation coefficients Q (Eq. 74 of Bucha et al., 2019)
%--------------------------------------------------------------------------
fprintf('Computing Q-terms (%s)\n',datestr(clock))

r_powers=zeros(1,max([pmax kmax])+1,'mp');
for e=0:max([pmax kmax])
    r_powers(e+1)=r^e;
end
Q=zeros(nmax+1,kmax+1,pmax,'mp');
for p_temp=1:pmax
    oo=0;
    for i=0:kmax
        if p_temp==1
            Q(:,oo+1,p_temp)=G(:,i+1);
        elseif p_temp==2
            Q(:,oo+1,p_temp)=-(i-1)*G(:,i+1)-r*G(:,i+2);
        else
            ap=(-1)^(p_temp-1)*fact(p_temp)*fact(p_temp-2);
            Q_temp=0;
            for s=1:p_temp-2
                aps=ap/fact(p_temp-s+1)/fact(p_temp-s-1)/fact(s);
                if i==0
                    Q_temp=Q_temp+aps*r_powers(p_temp-s+1)*G(:,p_temp-s+1);
                else
                    der_temp=0;
                    for k=0:i
                        ii=i-k;
                        if ii==0
                            der_temp=der_temp+binomial(i+1,k+1)*r_powers(p_temp-s+1)*G(:,p_temp-s+k+1);
                        else
                            if (i-k)>(p_temp-s)
                                continue
                            else
                                r_ampl=1;
                                for j=1:ii
                                    r_ampl=r_ampl*(p_temp-s-j+1);
                                end
                                r_ampl=r_ampl*r_powers(p_temp-s-ii+1);

                                der_temp=der_temp+binomial(i+1,k+1)*r_ampl*G(:,p_temp-s+k+1);
                            end
                        end
                    end
                    Q_temp=Q_temp+(aps*der_temp);
                end
            end
            Q(:,oo+1,p_temp)=Q_temp;
        end
        oo=oo+1;
    end
    Q(:,:,p_temp)=1/fact(p_temp+1).*Q(:,:,p_temp);
end
%--------------------------------------------------------------------------
