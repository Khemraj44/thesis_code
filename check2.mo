package ProjectHydroPower
//Model of generator connected to waterway and turbine
	model SimHydroPower
		//
		// Simulation experiment of hydro power system
		// Solution to project in course
		// FM1015 Modeling of Dynamic Systems
		// Fall semester 2016
		//
		// Bernt Lie
		// University College of Southeast Norway
		// September 26, 2016
		// modification of code by Liubomyr Vytvytskyi
		//
		ModHydroPower mhp; // Instantiating model of hydro power system
	//
		Real _v_f " Field voltage , V";
		//Real _T_m ;
		//parameter Real Tm0 = 85e4 ;//4e4;
		//parameter Real vf0 = 8500;//8500 ; //235;
		Real _u_v ;
		//Real _Wd_e ;
		Real _H_t ;
		Real _p_i1 , _p_m , _p_p2 , _p_tr2 , _p_d2 ;
		Real _w_g ;
		Real _h_s , _Vd_i , _Vd_s , _Vd_p ;
	equation
		_v_f = if time < 300 then 5000 else 4000;
		//_T_m = Tm0 ; //40000;
		mhp. v_f = _v_f ;
		//mhp. T_m = _T_m ;
		_u_v =if time < 300 then 0.5 else 0.47;
		//_Wd_e = if time < 1200 then 80 else 77;
		_H_t = 20; // if time < 1800 then 10 else 20;
	//
		mhp.u_v = _u_v ;
		//mhp. Wd_e = _Wd_e * 1e6;
		mhp.H_t = _H_t ;
		_p_i1 = mhp . p_i1 /mhp .p_a;
		_p_m = mhp.p_m/mhp. p_a ;
		_p_p2 = mhp . p_p2 /mhp .p_a;
		_p_tr2 = mhp. p_tr2 / mhp.p_a;
		_p_d2 = mhp . p_d2 /mhp .p_a;
		_w_g = 4* mhp. w_a /2/ mhp .PI;
		_h_s = mhp.h_s;
		_Vd_i = mhp . Vd_i ;
		_Vd_s = mhp . Vd_s ;
		_Vd_p = mhp . Vd_p ;
	end SimHydroPower ;
//
	model ModHydroPower
		//
		// Model of hydro power system
		// Solution to project in course
		// FM1015 Modeling of Dynamic Systems
		// Fall semester 2016
		// Constants
		constant Real g = 9.81 " Acceleration of gravity , m/s2 ";
		constant Real PI = 3.141592654 "pi , ";
		// General parameters
		parameter Real rho = 997 " Water density , kg/m3 ";
		parameter Real mu = 8.9e-4 " Viscosity of water , Pa.s";
		parameter Real epsilon = 1.5e-5 " Pipe roughness height , m";
		parameter Real H_r = 50 " Reservoir level above intake , m";
		parameter Real C_v = 6 " Guide vane capacity , m3/s";
		parameter Real p_a = 1.013e5 " Atmospheric pressure , Pa ";
		// Rated values
		parameter Real Vd_0 = 20 " Nominal flow rate , m3/s";
		//parameter Real w_a_0 = 50*2* PI /4 " Nominal angular velocity , rad/s";
		// Intake race
		parameter Real H_i = 25 " Vertical drop of intake race , m";
		parameter Real L_i = 6600 " Length of intake race , m";
		parameter Real D_i = 5.8 " Diameter of intake race , m";
		parameter Real A_i = PI*D_i ^2/4 " Cross section area of intake race , m2 ";
		parameter Real V_i = A_i* L_i " Volume of intake race , m3 ";
		parameter Real m_i = rho* V_i " Mass of intake race , kg ";
		parameter Real A_wi = PI* D_i*L_i " Wetting surface of intake race , m2 ";
		parameter Real p_i1 = p_a + rho*g*H_r " Intake race pressure ,Pa ";
		// Surge tank
		parameter Real H_s = 120 " Vertical drop of surge tank , m";
		parameter Real L_s = 140 " Length of surge tank , m";
		parameter Real D_s = 3.4 " Diameter of the surge tank , m";
		parameter Real A_s = PI*D_s ^2/4 " Cross section area of surge tank , m2 ";
		parameter Real p_s2 = p_a " Outlet pressure of surge tank , Pa";
		// Penstock
		parameter Real H_p = 420 " Vertical drop of penstock , m";
		parameter Real L_p = 600 " Length of penstock , m";
		parameter Real D_p = 3.3 " Diameter of penstock , m";
		parameter Real A_p = PI*D_p ^2/4 " Cross section area of penstock , m2 ";
		parameter Real V_p = A_p* L_p " Volume of penstock , m3 ";
		parameter Real m_p = rho* V_p " Mass of penstock , kg ";
		parameter Real A_wp = PI* D_p*L_p " Wetting surface of penstock , m2 ";
		// Discharge
		parameter Real H_d = 5 " Vertical drop of discharge race , m";
		parameter Real L_d = 600 " Length of discharge race , m";
		parameter Real D_d = 5.8 " Diameter of discharge race , m";
		parameter Real A_d = PI*D_d ^2/4 " Cross section area of discharge race , m2 ";
		parameter Real V_d = A_d* L_d " Volume of discharge race , m3 ";
		parameter Real m_d = rho* V_d " Mass of discharge race , kg ";
		parameter Real A_wd = PI* D_d*L_d " Wetting surface of discharge race , m2 ";
		// Aggregate
		parameter Real J = 3e6 " Moment of inertia of aggregate , kg .m2 ";
		parameter Real k_ba = 1e3 " Friction factor in the aggregate bearing box , W-s3/ rad3 ";
		parameter Real eta_e = 0.99 " Electricity generator efficiency ";
		// Turbine
		parameter Real eta_h = 0.90 " Turbine hydraulic efficiency ";
		//
		parameter Real Rf = 0.0715 " Field resistance , Ohm";
		parameter Real Lf = 576.92e-3 " Field inductance , H";
		parameter Real Mf = 32.653e-3 "Field - phase mutual inductance ,H";
		parameter Real Ra = 0.0031 " Phase a resistance , Ohm";
		parameter Real Rb = Ra " Phase b resistance , Ohm";
		parameter Real Rc = Ra " Phase c resistance , Ohm";
		parameter Real R = Ra " Common phase resistance , Ohm";
		parameter Real Lsl = 0.4129e-3 "armature leakage" ;
		parameter Real Ls = 3.2758e-3 + 1e-7 + Lsl " Static phase self inductance , H";
		parameter Real Lm = 0.0458e-3 " Magnitude phase self inductance , H";
		parameter Real Ms = 1.6379e-3 " Static phase - phase mutual inductance , H";
		parameter Real omega_0 = 2* PI *60 " Nominal grid electric angular velocity , rad/s";
		parameter Real V_L = 24e3 " Line voltage RMS , V";
		parameter Real V_phi = V_L / sqrt (3) " Phase voltage RMS , V";
		parameter Real Lg = 0.1 "Infinity bus inductance,H"; 
		// Load
		parameter Real delta0 = PI /6 " Load angle , rad";
		// Inductances
		parameter Real Ldd = Ls + Ms + 3* Lm /2;
		parameter Real Lqq = Ls + Ms - 3* Lm /2;
		parameter Real L00 = Ls - 2* Ms;
		//
		parameter Real Ldg = Ldd + Lg;
		parameter Real Lqg = Lqq + Lg;
		parameter Real Lff = Lf;
		//
		parameter Real Ldf = Mf* sqrt (3/2);
		parameter Real Lfd = Ldf;
		parameter Real Lqf = 0;
		parameter Real Lfq = Lqf;
		parameter Real L0f = 0;
		parameter Real Lf0 = L0f;
		// Initial state parameters
		parameter Real i_f0 = 1e3 " Field current initial value , A";
		parameter Real i_a0 = 1e4 " Phase a current initial value , A";
		parameter Real i_b0 = 1e4*sin ( -2* PI /3) " Phase b current initial value , A";
		parameter Real i_c0 = 1e4*sin ( -4* PI /3) " Phase c current initial value , A";
		parameter Real i_d0 = (cos (0) * i_a0 + cos ( -2* PI /3)* i_b0 + cos( -4* PI /3)* i_c0 )* sqrt (2/3)" Direct axis current initial value , A";
		parameter Real i_q0 = (sin (0) * i_a0 + sin ( -2* PI /3)* i_b0 + sin( -4* PI /3)* i_c0 )* sqrt (2/3)" Quadrature axis current initial value , A";
		parameter Real i_00 = ( i_a0 / sqrt (2) + i_b0 / sqrt (2) + i_c0 /sqrt (2))* sqrt (2/3)" Zero coordinate current initial value , A";
		parameter Real w_m0 = 2*PI*60;
		//
		parameter Real psi_d0 = Ldd * i_d0 + Ldf* i_f0 ;
		parameter Real psi_q0 = Lqq * i_q0 + Lqf* i_f0 ;
		parameter Real psi_00 = L00 * i_00 + L0f* i_f0 ;
		parameter Real psi_f0 = Lfd * i_d0 + Lfq* i_q0 + Lf0* i_00 + Lff*i_f0 ;
		// Variables
		Real m_s " Mass in surge tank , kg ";
		Real V_s " Volume in surge tank , m3 ";
		Real ell_s " Length of water string in surge tank , m";
		Real h_s " Level in surge tank , m";
		Real A_ws " Wetting surface of surge tank string , m2 ";
		//
		Real K3p_i , K3p_s , K3p_p , K3p_d ;// Volumetric kinetic energies
		Real md_i , md_s , md_p , md_d ; // mass flow rates
		Real Vd_i , Vd_s , Vd_p , Vd_d ; // Volumetric flow rates
		Real v_i , v_p , v_s , v_d ; // linear velocities
		Real N_Re_s , N_Re_i , N_Re_p , N_Re_d ; // Reynold numbers
		Real f_Di , f_Ds , f_Dp , f_Dd ; // Darcy friction coefficients
		Real F_fs , F_fi , F_fp , F_fd ; // Friction forces
		Real F_gs , F_gi , F_gp , F_gd ; // Gravity forces
		Real M_i , M_s , M_p , M_d ; // Linear momenta
		Real Md_s ; // Linear momentum flow into /out of surge tank
		Real p_m " Pressure in manifold , Pa ";
		Real p_i2 , p_s1 , p_p1 , p_p2 , p_tr2 , p_d1 , p_d2 ; // Pipe end pressures
		Real dp_tr " Pressure drop across turbine , Pa ";
		//Real K_a; // Aggregate kinetic energy
		// Real Kd_p , Kd_tr2 , Kd_d1 ; // Kinetic flow rates
		Real Wd_fa , Wd_ts ; // Work powers
		//Real w_m; // Angular velocity of aggregate
		// Variables
		Real i_f( start =i_f0 , fixed = true ) " Field current , A";
		Real psi_f ( start = psi_f0 ) " Field flux linkage , V.s";
		Real i_d( start =i_d0 , fixed = true ) " Direct axis current , A";
		Real psi_d ( start = psi_d0 ) " Direct axis flux linkage , V.s";
		Real i_q( start =i_q0 , fixed = true ) " Quadrature axis current , A";
		Real psi_q ( start = psi_q0 ) " Quadrature axis flux linkage , V.se";
		Real i_0( start =i_00 , fixed = true ) " Zero coordinate current , A";
		Real psi_0 ( start = psi_00 ) " Zero coordinate flux linkage , V.s";
		Real w_a (start = w_m0, fixed = true);
		Real delta(start = delta0, fixed = true);
		//
		Real omega_e " Electric angular velocity , rad/s";
		Real theta_e " Electric angle phase a, rad";
		Real theta_a " Electric angle phase a, rad";
		Real theta_b " Electric angle phase b, rad";
		Real theta_c " Electric angle phase c, rad";
		//
		Real v_ddd " Phase a voltage , V";
		Real v_q " Phase b voltage , V";
		Real v_0 " Phase c voltage , V";
		//
		Real i_a;
		Real i_b;
		Real i_c;
		//
		Real v_a;
		Real v_b;
		Real v_c;
		//
		Real v_a_g;
		Real v_b_g;
		Real v_c_g;
		// Torque
		Real T_e "Electromagnetic Torque, Nm"; 
		//
		Real i_d_ss ;
		Real i_q_ss ;
		Real den_i_q_ss ;
		//
		Real P " Power in the phases , W";
		Real Pr " Power loss in rotor , W";
		Real Ps " Power loss in stator , W";
		Real Pa " Apparent power , VA";
		Real eta_P " Power factor ,_";
		Real P_dq0 " Power in the phases in dq0 coordinates , W";
		Real Ps_dq0 " Power loss in stator in dq0 coordinates , W";
		Real P_ss " Steady state power , W";
		Real Pr_ss " Steady state power loss in rotor , W";
		Real Ps_ss " Steady state power loss in stator , W";
		Real Pa_ss " Steady state apparent power , VA";
		Real eta_P_ss " Steady state power factor , -";
		Real I_phi " Steady state phase current amplitude , A";
		//Real P_m "mechanical power";
		input Real v_f " Field voltage , V";
		//input Real T_m ;
		// Inputs
		input Real u_v " Input guide vane signal , -";
		Real Wd_e;
		//input Real Wd_e " Electrical power usage , W";
		input Real H_t " Tail water level , m";
		output Real y_P;
	initial equation
		Vd_i = Vd_0 ;
		Vd_p = Vd_0 ;
		//w_m = w_m0 ;
		h_s = H_r + H_i;
	equation
		// Surge tank liquid string
		m_s = rho *V_s;
		V_s = A_s * ell_s ;
		h_s = ell_s *H_s /L_s;
		A_ws = PI* D_s* ell_s " Wetting surface of surge tank string ,m2 ";
		// Discharge counter pressure
		p_d2 = p_a + rho * g * H_t ;
		// Mass flow rates
		md_i = rho* Vd_i ;
		md_s = rho* Vd_s ;
		md_p = rho* Vd_p ;
		md_d = rho* Vd_d ;
		// Linear velocities
		Vd_i = A_i*v_i;
		Vd_s = A_i*v_s;
		Vd_p = A_p*v_p;
		Vd_d = A_d*v_d;
		// Manifold
		p_i2 = p_m;
		p_s1 = p_m;
		p_p1 = p_m;
		Vd_i = Vd_s + Vd_p ;
		// Penstock - Discharge string
		Vd_d = Vd_p ;
		// Momentums
		M_i = m_i *v_i;
		M_s = m_s *v_s;
		M_p = m_p *v_p;
		M_d = m_d *v_d;
		// Momentum flow rate into /out of surge tank
		Md_s = rho* Vd_s * abs( Vd_s )/A_s;
		// Gravitational forces
		F_gi = m_i*g;
		F_gs = m_s*g;
		F_gp = m_p*g;
		F_gd = m_d*g;
		// Reynold 's numbers
		N_Re_i = rho* abs(v_i)* D_i /mu;
		N_Re_s = rho* abs(v_s)* D_s /mu;
		N_Re_p = rho* abs(v_p)* D_p /mu;
		N_Re_d = rho* abs(v_d)* D_d /mu;
		// Darcy friction coefficients
		f_Di = fDarcy (N_Re_i , D_i , epsilon );
		f_Ds = fDarcy (N_Re_s , D_s , epsilon );
		f_Dp = fDarcy (N_Re_p , D_p , epsilon );
		f_Dd = fDarcy (N_Re_d , D_d , epsilon );
	// Volumetric kinetic energies
		K3p_i = rho *v_i*abs(v_i) /2;
		K3p_s = rho *v_s*abs(v_s) /2;
		K3p_p = rho *v_p*abs(v_p) /2;
		K3p_d = rho *v_d*abs(v_d) /2;
		// Friction forces
		F_fi = K3p_i * A_wi * f_Di /4;
		F_fs = K3p_s * A_ws * f_Ds /4;
		F_fp = K3p_p * A_wp * f_Dp /4;
		F_fd = K3p_d * A_wd * f_Dd /4;
		// Turbine
		Vd_p = C_v*u_v* sqrt ( dp_tr /p_a);
		dp_tr = p_p2 - p_tr2 ;
		Wd_ts = eta_h * dp_tr * Vd_p ;
		// Aggregate
		//K_a = J_a *w_m ^2/2;
		Wd_e = eta_e * P_dq0 ;
		Wd_fa = k_ba *w_a ^2/2;
		p_d1 = p_tr2 ;
		// Balance laws
		der(M_i) = ( p_i1 - p_i2 )*A_i + F_gi *H_i/ L_i - F_fi ;
		der(m_s) = md_s ;
		der(M_s) = Md_s + ( p_s1 - p_s2 )*A_s - F_gs *H_s/L_s - F_fs ;
		der(M_p) = (p_p1 - p_p2 )*A_p + F_gp *H_p/L_p - F_fp ;
		der(M_d) = (p_d1 - p_d2 )*A_d + F_gd *H_d/L_d - F_fd ;
		der(J *w_a ^2/2) = Wd_ts - Wd_fa - P_dq0 ;
		//
		// Generator..........................................................................................................
		omega_e = omega_0 ;
		theta_e = omega_e * time ;
		theta_a = theta_e + delta ;
		theta_b = theta_a - 2* PI /3;
		theta_c = theta_b - 2* PI /3;
		// Setting up voltages
		v_ddd = ( cos( theta_a )* v_a + cos( theta_b )*v_b + cos ( theta_c )*v_c)* sqrt (2/3) ;
		v_q = ( sin( theta_a )* v_a + sin( theta_b )*v_b + sin ( theta_c )*v_c)* sqrt (2/3) ;
		v_0 = ( v_a + v_b + v_c)/ sqrt (3);
		// Kirchhoff ’s voltage laws
		v_f = Rf* i_f + der( psi_f );
		v_ddd = -R* i_d - der( psi_d ) - omega_e * psi_q ;
		v_q = -R* i_q - der( psi_q ) + omega_e * psi_d ;
		v_0 = -R* i_0 - der( psi_0 );
		// Flux linkage vs current
		psi_d = Ldd *i_d + Ldf *i_f ;
		psi_q = Lqq *i_q + Lqf *i_f ;
		psi_0 = L00 *i_0 + L0f *i_f ;
		psi_f = Lfd *i_d + Lfq *i_q + Lf0 *i_0 + Lff*i_f;
		// Phase currents
		i_a = ( cos( theta_a )* i_d + sin( theta_a )*i_q + i_0 / sqrt (2) )*sqrt (2/3) ;
		i_b = ( cos( theta_b )* i_d + sin( theta_b )*i_q + i_0 / sqrt (2) )*sqrt (2/3) ;
		i_c = ( cos( theta_c )* i_d + sin( theta_c )*i_q + i_0 / sqrt (2) )*sqrt (2/3) ;
		// Phase grid voltages
		v_a_g = V_phi * sqrt (2)* sin( theta_e );
		v_b_g = V_phi * sqrt (2)* sin( theta_e - 2 * PI /3) ;
		v_c_g = V_phi * sqrt (2)* sin( theta_e - 4 * PI /3) ;
		//phase voltages
		v_a = Lg * der(i_a) + v_a_g ;
		v_b = Lg * der(i_b) + v_b_g ;
		v_c = Lg * der(i_c) + v_c_g ;
		//Torque
		T_e = psi_d*i_q - psi_q*i_d;
		//Swing equation
		//Pm0 = T_m*w_m0;
		//der(w_m) = 1/J * (T_m - T_e ) ;//replace with 1/J * (T_m - T_e - kf*w_m)  //- 8488.26*(w_m-omega_e)
		der(delta) = w_a - omega_e;
		// Power
		//P_m = T_m*w_m;
		P = v_a *i_a + v_b *i_b + v_c *i_c ;
		P_dq0 = v_ddd *i_d + v_q *i_q + v_0 *i_0;
		Pr = Rf*i_f ^2;
		Ps = Ra*i_a ^2 + Rb*i_b ^2 + Rc*i_c ^2;
		Ps_dq0 = R*( i_d ^2 + i_q ^2 + i_0 ^2);
		Pa = V_L* sqrt (i_d ^2 + i_q ^2);
		eta_P = P_dq0 /Pa;
		//
		den_i_q_ss = R^2 + omega_e ^2* Ldg *Lqg ;
		i_d_ss = -( sqrt (3/2) * omega_e ^2* Lqg*Mf* i_f - (R* sin( delta )+omega_e *Lqq*cos( delta ))* V_L)/ den_i_q_ss ;
		i_q_ss = ( sqrt (3/2) * omega_e *Mf*R* i_f - (R* cos( delta )-omega_e *Ldg* sin( delta ))* V_L)/ den_i_q_ss ;
		//
		P_ss = V_L *( i_q_ss *cos( delta ) - i_d_ss * sin( delta ));
		Pr_ss = v_f ^2/ Rf;
		Ps_ss = R*( i_d_ss ^2 + i_q_ss ^2) ;
		Pa_ss = V_L * sqrt ( i_d_ss ^2 + i_q_ss ^2);
		eta_P_ss = P_ss / Pa_ss ;
		//
		I_phi = sqrt (2/3) * sqrt ( i_d_ss ^2 + i_q_ss ^2);
		//
		y_P = P;	
	end ModHydroPower ;
		//
		function fDarcy " Darcy friction factor "
		// Function for computing Darcy 's friction factor
		// author : Bernt Lie
		// University College of Southeast Norway
		// September 26, 2016
		//
		// Function input arguments
		input Real N_Re " Reynold number , -";
		input Real D " Pipe diameter , m";
		input Real epsilon " Pipe roughness height , m";
		// Function output ( response ) value
		output Real fD " Darcy friction factor , -";
		// Local ( protected ) quantities
		protected
		Real arg;
		// Algorithm for computing specific enthalpy
	algorithm
		arg := epsilon /3.7/ D + 5.74/(N_Re + 1e-3)^0.9;
		if arg <= 0 then
			fD :=0;
		else
			fD := 1/(2* log10 ( arg))^2;		
		end if;
	end fDarcy ;
//
end ProjectHydroPower ;