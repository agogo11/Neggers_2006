PRO yzreproduce_fig4_neggers_2006
  ;+
  ;	NAME:
  ;			yzReproduce_Fig4_Neggers_2006
  ; PURPOSE:
  ;     This is an IDL code to reproduce the Figure 4 in Neggers et al. (2006)
  ;
  ;     Neggers, R., Stevens, B., & Neelin, J. D. (2006).
  ;     A simple equilibrium model for shallow-cumulus-topped mixed layers.
  ;     Theoretical and Computational Fluid Dynamics, 20(5-6), 305-322.
  ;
  ; EXPLANATION:
  ;
  ;     The plot requires the coyote IDL package
  ;
  ;	CALLING SEQUENCE:
  ;
  ;	INPUTS:
  ;
  ;	OPTIONAL INPUT KEYWORDS:
  ;
  ;	OUTPUTS:
  ;
  ;	OPTIONAL OUTPUT KEYWORDS:
  ;
  ;	RESTRICTIONS:
  ;	    -The output are a little bit different from Neggers et al. (2006).
  ;	    -In particular, the calculated a_c in equiligrium doubles that from Neggers et al. (2006)
  ;	    -Reasons have not been found.
  ;
  ; AUTHOR:
  ; 		Youtong Zheng, Post-doc Associate
  ;			Department of Atmospheric and Oceanic Science &
  ;			Earth System Science Interdisciplinary Center (ESSIC)
  ;			University of Maryland
  ;			5825 University Research Court, Suite 4001
  ;			College Park, MD 20740-3823
  ;			Phone: 1-301-526-5086
  ;			Email: ytzheng@umd.edu
  ;
  ;	MODIFICATION HISTORY:
  ;			Written by: Youtong Zheng, Aug 2, 2018
  ;-

  Path = 'G:\MyIDLpackage\MLM_Neggers2006\'

  ;constants
  g=9.80665d ; gravity const ms-2
  BOLTZMAN = 1.380658e-23;            % BOLTZMAN,  Boltzman constant (J/K)
  AVOGADRO = .602214199e24;           % AVOGADRO, Avogadro constant (1/mol)
  MD       = 28.9644e-3;              % MD, molar mass dry air (kg/mol)
  MV  = 18.0153e-3;                   % MD, molae mass water vapor (kg/mol)
  r_v = (AVOGADRO)*(BOLTZMAN)/(MV); % r_v, gas constant for water vapor (J/(kg-K))
  r_d = (AVOGADRO)*(BOLTZMAN)/(MD); % r_d, gas constant for dry air (J/(kg-K))
  cp  = 7./2d*(r_d);                   % cp, specific heat of air (J/(kg-K))
  lv  = 2.5008e6;                     % lv, latent heat of vaporization (J/kg)
  rho = 1.d ;kg/m3

  ;Set input parameters
  V = 8.75; surface wind speed m/s
  SST = 300.4; sea surface temperature K
  D = 4.3e-6; large-scale divergence rate s-1
  qt_plus = 4.0; free-atmosphere water vapor mixing ratio g/kg
  qt_plus = qt_plus/1000. ; g/kg -> kg/kg
  thetal_plus = 308.; free-troposphere potential temperature K

  F_rad = -2.; mixed-layer radiative cooling Kday-1
  F_rad = F_rad/(24.*60.*60d) ; Ks-1

  F_adv_thetal = 0d; potential temperature horizontal advection Kday-1

  F_adv_qt = -1.2d; moisture advection gkg-1day-1
  F_adv_qt = (F_adv_qt/1000d)/(24.*60.*60d) ; kgkg-1s-1

  ps = 1015d; surface presure hPa
  p_ref = 1000.; reference pressure

  ;Empirical coefficients
  C_theta_s = 0.0012d
  C_q_s = 0.0012d

  C_theta_c = 0.03d
  C_q_c = 0.1d

  deltaZ = 100d ; m

  ;Integration time step
  deltat = 60d ;s

  ;initial condition
  thetal_s0 = theta(ps,SST - 273.15d)
  Ts0 = thetal_s0/((1000./ps)^0.28571429)
  thetal_m0 = thetal_s0 - 2d
  T20 = thetal_m0/((1000./ps)^0.28571429)
  qt_m0 = qt_plus

  qt_s0 = mixr_sat(SST, ps)
  qt_s0 = 0.981*qt_s0/1000d

  h0 = [100d,400d,700d,1000d]; m

  thetal_s = thetal_s0
  qt_s = qt_s0

  M0 = 0
  ;three prognostic variables
  qt_m = qt_m0
  thetal_m = thetal_m0

  ;output parameters
  nt = 167*60. + 1
  t_plot = (1./60)*FINDGEN(nt)
  h_plot = FLTARR(4, nt)
  a_c_plot = FLTARR(4, nt)
  M_plot = FLTARR(4, nt)
  E_plot = FLTARR(4, nt)
  qt_m_plot = FLTARR(4, nt)
  LHF_plot = FLTARR(4, nt)
  sigma_q_plot = FLTARR(4, nt)
  zb_h_plot = FLTARR(4, nt)

  dqt_m_residual = FLTARR(4, nt)
  dqt_m_LHF = FLTARR(4, nt)
  dqt_m_E = FLTARR(4, nt)
  dqt_m_adv = FLTARR(4, nt)

  FOR icases = 0, 1 DO BEGIN
    SST = 300.4
    thetal_s = theta(ps,SST - 273.15d)
    thetal_m = thetal_s0 - 2d

    qt_m = qt_plus

    qt_s = mixr_sat(SST, ps)
    qt_s = 0.981*qt_s/1000d

    h = 100
    M = 0

    D = 4.3e-6

    FOR t = 0L, 601200L, 60L DO BEGIN

      IF t EQ 117L*60*60 AND icases EQ 1 THEN BEGIN
        PRINT, 'perturb 1 kg/kg'
        qt_m = qt_m + 1./1000.
      ENDIF

      ;large-scale vertical velocity
      w = -D*h

      ;surface heat fluxes
      flux_thetal_s = C_theta_s*V*(thetal_s - thetal_m)
      flux_qt_s = C_q_s*V*(qt_s - qt_m)

      flux_thetav_s = flux_thetal_s*(1. + 0.61*((qt_m + qt_s)/2.)) + 0.61*((thetal_s + thetal_m)/2.)*flux_qt_s

      ;inversion
      delta_qt = C_q_c*(qt_plus - qt_m)
      delta_thetal = C_theta_c*(thetal_plus - thetal_m)

      delta_thetav = delta_thetal + 0.61*(qt_m*delta_thetal + thetal_m*delta_qt + delta_qt*delta_thetal)

      ;entrainment
      E = 0.2*flux_thetav_s/delta_thetav

      ;cloud-base mass flux
      thetav_m = thetal_m*(1. - 1.61*qt_m)

      IF t EQ 0 THEN thetav_m0 = thetav_m

      ww = ((g*h/thetav_m)*flux_thetav_s)^(1./3.)
      w_c = ww

      sigma_q = SQRT(-(flux_qt_s/w_c)*delta_qt*(h/deltaz))
      ;sigma_q = abs(-(flux_qt_s/w_c)*delta_qt*(h/deltaz))

      ;surface air temp
      T2 = thetal_m/((p_ref/ps)^0.28571429)
      qt_s_sat = mixr_sat(T2, ps)/1000. ;kg/kg

      ;mixed-layer top temp
      p_h = ps/EXP((g*h)/(r_d*thetal_m))
      T_h = thetal_m/((ps/p_h)^(r_d/cp))

      ;moisture deficit
      qt_b_sat = mixr_sat(T_h + delta_thetal, p_h)/1000. ;kg/kg
      qt_b = qt_m + delta_qt
      defi = qt_b - qt_b_sat

      a_c = 0.5d + 0.36d*ATAN(1.55d*defi/sigma_q)

      M = a_c*w_c

      ;other related parameters
      ;zb - h
      gradient_qt_sat = (qt_s_sat - qt_b_sat)/h
      zb_h = -defi/gradient_qt_sat

      IF a_c LE 0. THEN M = 0

      ;prognostic variables
      h_new = h + deltat*(E + w - M)
      qt_m_new = qt_m + (deltat/h)*(flux_qt_s + E*delta_qt + h*F_adv_qt)
      thetal_m_new = thetal_m + (deltat/h)*(flux_thetal_s + E*delta_thetal + h*F_adv_thetal + h*F_rad)

      ;others
      LHF = rho*Lv*flux_qt_s
      RH_s = mixr2rh(1000.*qt_m,ps,SST - 2)
      LCL = romplcl(ps*100,SST - 2,RH_s/100.)

      it = t/60.

      h_plot[icases,it] = h
      a_c_plot[icases,it] = a_c
      M_plot[icases,it] = M
      E_plot[icases,it] = E
      qt_m_plot[icases,it] = qt_m
      LHF_plot[icases,it] = LHF
      sigma_q_plot[icases,it] = sigma_q
      zb_h_plot[icases,it] = LCL - h
      zb_h_plot[icases,it] = zb_h

      h = h_new
      qt_m = qt_m_new
      thetal_m = thetal_m_new

      ;PRINT, h_new, E, a_c, M

      ;PRINT, ' '
    ENDFOR
  ENDFOR

  PRINT, ' '

  cgps_open, Path + 'Fig_Neggers_2006_Fig4_new.eps',/CMYK,$
    /nomatch,Font = 1,/Times, xsize =12/2.54, ysize = 10/2.54, fontsize =14

  pos = cglayout([2,4], OXMargin=[4.,1.], OYMargin=[3.,1.],Ygap = 2.5, Xgap = 4.)

  ;Fig 1
  cgplot,t_plot, h_plot,/noerase,/nodata, $
    xtitle = '', ytitle = 'h [m]',$
    xrange = [0.,167.],yrange = [0.,2000],$
    charsize = 0.6, POSITION = pos[*,0]

  cgoplot,t_plot, h_plot[0,*], linestyle = 0, thick = 1.
  cgoplot,t_plot, h_plot[1,*], linestyle = 1, thick = 1.
  cgoplot,t_plot, h_plot[2,*], linestyle = 2, thick = 1.
  cgoplot,t_plot, h_plot[3,*], linestyle = 3, thick = 1.

  ;Fig 2
  cgplot,t_plot, h_plot,/noerase,/nodata, $
    xtitle = '', ytitle = 'qt_m [g kg-1]',$
    xrange = [0.,167.],yrange = [0.,20],$
    charsize = 0.6, POSITION = pos[*,1]

  cgoplot,t_plot, 1000.*qt_m_plot[0,*], linestyle = 0, thick = 1.
  cgoplot,t_plot, 1000.*qt_m_plot[1,*], linestyle = 1, thick = 1.
  cgoplot,t_plot, 1000.*qt_m_plot[2,*], linestyle = 2, thick = 1.
  cgoplot,t_plot, 1000.*qt_m_plot[3,*], linestyle = 3, thick = 1.

  ;Fig 3
  cgplot,t_plot, h_plot,/noerase,/nodata, $
    xtitle = '', ytitle = 'a_c [%]',$
    xrange = [0.,167.],yrange = [0.,20],$
    charsize = 0.6, POSITION = pos[*,2]

  cgoplot,t_plot, 100.*a_c_plot[0,*], linestyle = 0, thick = 1.
  cgoplot,t_plot, 100.*a_c_plot[1,*], linestyle = 1, thick = 1.
  cgoplot,t_plot, 100.*a_c_plot[2,*], linestyle = 2, thick = 1.
  cgoplot,t_plot, 100.*a_c_plot[3,*], linestyle = 3, thick = 1.

  ;Fig 4
  cgplot,t_plot, h_plot,/noerase,/nodata, $
    xtitle = '', ytitle = 'LHF [Wm-2]',$
    xrange = [0.,167.],yrange = [0.,500],$
    charsize = 0.6, POSITION = pos[*,3]

  cgoplot,t_plot, LHF_plot[0,*], linestyle = 0, thick = 1.
  cgoplot,t_plot, LHF_plot[1,*], linestyle = 1, thick = 1.
  cgoplot,t_plot, LHF_plot[2,*], linestyle = 2, thick = 1.
  cgoplot,t_plot, LHF_plot[3,*], linestyle = 3, thick = 1.

  ;Fig 5
  cgplot,t_plot, h_plot,/noerase,/nodata, $
    xtitle = '', ytitle = 'M [ms-1]',$
    xrange = [0.,167.],yrange = [0.,0.1],$
    charsize = 0.6, POSITION = pos[*,4]

  cgoplot,t_plot, M_plot[0,*], linestyle = 0, thick = 1.
  cgoplot,t_plot, M_plot[1,*], linestyle = 1, thick = 1.
  cgoplot,t_plot, M_plot[2,*], linestyle = 2, thick = 1.
  cgoplot,t_plot, M_plot[3,*], linestyle = 3, thick = 1.

  ;Fig 6
  cgplot,t_plot, h_plot,/noerase,/nodata, $
    xtitle = '', ytitle = 'sigma_q [g kg-1]',$
    xrange = [0.,167.],yrange = [0.,1.5],$
    charsize = 0.6, POSITION = pos[*,5]

  cgoplot,t_plot, 1000.*sigma_q_plot[0,*], linestyle = 0, thick = 1.
  cgoplot,t_plot, 1000.*sigma_q_plot[1,*], linestyle = 1, thick = 1.
  cgoplot,t_plot, 1000.*sigma_q_plot[2,*], linestyle = 2, thick = 1.
  cgoplot,t_plot, 1000.*sigma_q_plot[3,*], linestyle = 3, thick = 1.

  ;Fig 7
  cgplot,t_plot, h_plot,/noerase,/nodata, $
    xtitle = 'Time [hr]', ytitle = 'E [ms-1]',$
    xrange = [0.,167.],yrange = [0.,0.1],$
    charsize = 0.6, POSITION = pos[*,6]

  cgoplot,t_plot, E_plot[0,*], linestyle = 0, thick = 1.
  cgoplot,t_plot, E_plot[1,*], linestyle = 1, thick = 1.
  cgoplot,t_plot, E_plot[2,*], linestyle = 2, thick = 1.
  cgoplot,t_plot, E_plot[3,*], linestyle = 3, thick = 1.

  ;Fig 8
  cgplot,t_plot, h_plot,/noerase,/nodata, $
    xtitle = 'Time [hr]', ytitle = 'zb - h [m]',$
    xrange = [0.,167.],yrange = [0.,1500],$
    charsize = 0.6, POSITION = pos[*,7]

  cgoplot,t_plot, zb_h_plot[0,*], linestyle = 0, thick = 1.
  cgoplot,t_plot, zb_h_plot[1,*], linestyle = 1, thick = 1.
  cgoplot,t_plot, zb_h_plot[2,*], linestyle = 2, thick = 1.
  cgoplot,t_plot, zb_h_plot[3,*], linestyle = 3, thick = 1.

  cgps_close, /png, density = 1000

  PRINT,''

END