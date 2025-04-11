#pragma once
// contactmodellinearpbond.h

#include "contactmodel/src/contactmodelmechanical.h"

#ifdef macroelement_LIB
#  define macroelement_EXPORT EXPORT_TAG
#elif defined(NO_MODEL_IMPORT)
#  define macroelement_EXPORT
#else
#  define macroelement_EXPORT IMPORT_TAG
#endif

namespace cmodelsxd {
	using namespace itasca;

	class ContactModelmacroelement : public ContactModelMechanical {
	public:
		macroelement_EXPORT ContactModelmacroelement();
		macroelement_EXPORT virtual ~ContactModelmacroelement();

		virtual QString  getName() const { return "macroelement"; }
		virtual void     setIndex(int i) { index_ = i; }
		virtual int      getIndex() const { return index_; }
		virtual void                     copy(const ContactModel *c);
		virtual void                     archive(ArchiveStream &);


		enum PropertyKeys {
			kwLinKn = 1
			, kwLinKs
			, kwLinFric
			, kwLinF
			, kwLinS
			, kwLinMode
			, kwRGap
			, kwEmod
			, kwKRatio
			, kwDpNRatio
			, kwDpSRatio
			, kwDpMode
			, kwDpF
			, kwPbState
			, kwPbRMul
			, kwPbKn
			, kwPbKs
			, kwPbMcf
			, kwPbTStrength
			, kwPbSStrength
			, kwPbCoh
			, kwPbFa
			, kwPbSig
			, kwPbTau
			, kwPbF
			, kwPbM
			, kwPbRadius
			, kwPbEmod
			, kwPbKRatio
			, kwUserArea
			, kwEA
			, kwEI
			, kwD
			, kwRefL
			, kwCurrL
			, kwNZero
			, kwNult
			, kwNtens		
			, kwMff
			, kwBetaf
			, kwKF
			, kwKM
			, kwMfg
			, kwAlphag
			, kwBetag
			, kwWpl1
			, kwWpl2
			, kwDalpha
			, kwRatA
			, kwRatB
			, kwNf
			, kwUP
			, kwOmegaP
			, kwForce
			, kwMoment
			, kwIsPlastic
			, kwMomentSign
			, kwUSlide
			, kwWP
			, kwForceDamaged
			, kwMomentDamaged

		};

		virtual QString  getProperties() const {
			return "kn"
				",ks"
				",fric"
				",lin_force"
				",lin_slip"
				",lin_mode"
				",rgap"
				",emod"
				",kratio"
				",dp_nratio"
				",dp_sratio"
				",dp_mode"
				",dp_force"
				",pb_state"
				",pb_rmul"
				",pb_kn"
				",pb_ks"
				",pb_mcf"
				",pb_ten"
				",pb_shear"
				",pb_coh"
				",pb_fa"
				",pb_sigma"
				",pb_tau"
				",pb_force"
				",pb_M"
				",pb_radius"
				",pb_emod"
				",pb_kratio"
				",user_area"
				",me_ea"
				",me_ei"
				",me_d"
				",me_ref_l"
				",me_curr_l"
				",me_nzero"
				",me_nult"
				",me_ntens"
				",me_mff"
				",me_betaf"
				",me_kf"
				",me_km"
				",me_mfg"
				",me_alphag"
				",me_betag"
				",me_wpl1"
				",me_wpl2"
				",me_dalpha"
				",me_rata"
				",me_ratb"
				",me_nf"
				",me_up"
				",me_omegap"
				",me_force"
				",me_moment"
				",me_isplastic"
				",me_moment_sign"
				",me_uslide"
				",me_wp"
				",me_force_damaged"
				",me_moment_damaged"
				;

		}

		enum EnergyKeys { kwEStrain = 1, kwESlip, kwEDashpot, kwEPbStrain, kwMeEStrain, kwMeEPlastic, kwMeEBStrain, kwMeEBPlastic};
		virtual QString  getEnergies() const { return "energy-strain,energy-slip,energy-dashpot,energy-pbstrain,energy-mestrain,energy-meplastic,energy-mebstrain,energy-mebplastic"; }
		virtual double   getEnergy(uint i) const;  // Base 1
		virtual bool     getEnergyAccumulate(uint i) const; // Base 1
		virtual void     setEnergy(uint i, const double &d); // Base 1
		virtual void     activateEnergy() { if (energies_) return; energies_ = NEWC(Energies()); }
		virtual bool     getEnergyActivated() const { return (energies_ != 0); }


		enum FishCallEvents { fActivated = 0, fBondBreak, fSlipChange };
		virtual QString  getFishCallEvents() const { return "contact_activated,bond_break,slip_change"; }
		virtual QVariant getProperty(uint i, const IContact *con = 0) const;
		virtual bool     getPropertyGlobal(uint i) const;
		virtual bool     setProperty(uint i, const QVariant &v, IContact *con = 0);
		virtual bool     getPropertyReadOnly(uint i) const;

		virtual bool     supportsInheritance(uint i) const;
		virtual bool     getInheritance(uint i) const { assert(i < 32); quint32 mask = to<quint32>(1 << i);  return (inheritanceField_ & mask) ? true : false; }
		virtual void     setInheritance(uint i, bool b) { assert(i < 32); quint32 mask = to<quint32>(1 << i);  if (b) inheritanceField_ |= mask;  else inheritanceField_ &= ~mask; }

		enum MethodKeys {
			kwDeformability = 1
			, kwPbDeformability
			, kwPbBond
			, kwPbUnbond
			, kwArea
		};

		virtual QString  getMethods() const {
			return "deformability"
				",pb_deformability"
				",bond"
				",unbond"
				",area";
		}

		virtual QString  getMethodArguments(uint i) const;

		virtual bool     setMethod(uint i, const QVector<QVariant> &vl, IContact *con = 0); // Base 1 - returns true if timestep contributions need to be updated

		virtual uint     getMinorVersion() const;

		virtual bool    validate(ContactModelMechanicalState *state, const double &timestep);
		virtual bool    endPropertyUpdated(const QString &name, const IContactMechanical *c);
		virtual bool    forceDisplacementLaw(ContactModelMechanicalState *state, const double &timestep);
		virtual DVect2  getEffectiveTranslationalStiffness() const { DVect2 ret = effectiveTranslationalStiffness_; if (pbProps_) ret += pbProps_->pbTransStiff_; return ret; }
		virtual DAVect  getEffectiveRotationalStiffness() const { if (!pbProps_) return DAVect(0.0); return pbProps_->pbAngStiff_; }

		struct systemStatus {
			//	systemStatus() : u_p_(0.0), theta_p_(0.0), force_(0.0), M_(0.0), alpha_(0.0) {}
			double u_p_;
			double omega_p_;
			double force_;
			double moment_;
			double u_slide_;
			double wp_;
			double nf_;
			double gap_;
			int moment_sign_ = 1;
		};

		struct epIncrements {
			//	epIncrements() : force_dot_(0.0), M_dot_(0.0), theta_p_dot_(0.0), u_p_dot_(0.0), alpha_dot_(0.0) {}
			double force_dot_;
			double moment_dot_;
			double u_p_dot_;
			double omega_p_dot_;
			double nf_dot_;
			double wp_;
		};

		struct dfIncrements {
			//	dfIncrements() : theta_dot_(0.0), u_dot_(0.0) {}
			double u_dot_;
			double omega_dot_;
		};

		// my functions
		virtual systemStatus stressStateEvolution(systemStatus current_state,DVect trans, DAVect ang, double gap, ContactModelMechanicalState* state);	// deals with forces/Ms in the DVect PFC format and evaluate the plastic status with yieldf
		virtual epIncrements evaluatePlasticity(systemStatus mystatus, dfIncrements dfinc);					// calculates the plastic increments given system status and strains using the constitutive model
		virtual systemStatus plasticIncrements(systemStatus curr_state, dfIncrements dfinc);				// substepping scheme for the plastic equation integration using RK23
		virtual double yieldf(double u_p,double theta_p, double force, double moment);						// evaluate the position of the current stress state in relation to the yield function
		virtual systemStatus updateStatus(systemStatus myStatus, epIncrements myIncrements);				// add the plastic increments to the system status
		virtual epIncrements epIncProduct(epIncrements myIncrements, double myMult);						// multiplies the plastic increments by a double, used for Runge-Kutta
		virtual epIncrements epIncSum(epIncrements myIncrements1, epIncrements myIncrements2);				// sums the plastic increments, used for Runge-Kutta
		virtual double calculateElasticModulus(double f);													// get the elastic stiffness given the current force


		virtual bool thermalCoupling(ContactModelMechanicalState *, ContactModelThermalState *, IContactThermal *, const double &);

		virtual ContactModelmacroelement *clone() const { return NEWC(ContactModelmacroelement()); }
		virtual double  getActivityDistance() const { return rgap_; }
		virtual bool    isOKToDelete() const { return !isBonded(); }
		virtual void    resetForcesAndMoments() { lin_F(DVect(0.0)); dp_F(DVect(0.0)); pbF(DVect(0.0)); pbM(DAVect(0.0)); me_force(0.0); 
		if (energies_) { energies_->estrain_ = 0.0; energies_->me_estrain_ = 0.0; energies_->me_eplastic_ = 0.0; energies_->me_ebstrain_ = 0.0; energies_->me_ebplastic_ = 0.0; energies_->epbstrain_ = 0.0; } }
		virtual void    setForce(const DVect &v, IContact *c);
		virtual void	setArea(const double &d) { userArea_ = d; }
		virtual double	getArea() const { return userArea_; } //added in PFC7

		virtual bool     checkActivity(const double &gap) { return (gap <= rgap_ || isBonded()); }

		virtual bool     isSliding() const { return lin_S_; }
		virtual bool     isBonded() const { return pbProps_ ? (pbProps_->pb_state_ == 3) : false; } // still using the parallel bond methods
		virtual void     propagateStateInformation(IContactModelMechanical* oldCm, const CAxes &oldSystem = CAxes(), const CAxes &newSystem = CAxes());
		virtual void     setNonForcePropsFrom(IContactModel *oldCM);

		const double & kn() const { return kn_; }
		void           kn(const double &d) { kn_ = d; }
		const double & ks() const { return ks_; }
		void           ks(const double &d) { ks_ = d; }
		const double & fric() const { return fric_; }
		void           fric(const double &d) { fric_ = d; }
		const DVect &  lin_F() const { return lin_F_; }
		void           lin_F(const DVect &f) { lin_F_ = f; }
		bool           lin_S() const { return lin_S_; }
		void           lin_S(bool b) { lin_S_ = b; }
		uint           lin_mode() const { return lin_mode_; }
		void           lin_mode(uint i) { lin_mode_ = i; }
		const double & rgap() const { return rgap_; }
		void           rgap(const double &d) { rgap_ = d; }

		// elastic model parameters
		const double & me_ea() const { return me_ea_; }
		void		   me_ea(const double &d) { me_ea_ = d;}
		const double& me_ei() const { return me_ei_; }
		void		   me_ei(const double& d) { me_ei_ = d; }

		const double& me_d() const { return me_d_; }
		void		   me_d(const double& d) { me_d_ = d; }
		const double& me_ref_l() const { return me_ref_l_; }
		void		   me_ref_l(const double& d) { me_ref_l_ = d; }
		const double& me_curr_l() const { return me_curr_l_; }
		void		   me_curr_l(const double& d) { me_curr_l_ = d; }

		// yield surface parameters
		const double& me_mff() const { return me_mff_; }
		void		   me_mff(const double& d) { me_mff_ = d; }
		const double& me_betaf() const { return me_betaf_; }
		void		   me_betaf(const double& d) { me_betaf_ = d; }
		const double & me_nzero() const { return me_nzero_; }
		void		   me_nzero(const double &d) { me_nzero_ = d; }
		const double & me_nult() const { return me_nult_; }
		void		   me_nult(const double &d) { me_nult_ = d; }
		const double & me_ntens() const { return me_ntens_; }
		void		   me_ntens(const double &d) { me_ntens_ = d; }

		// hardening function parameters
		const double & me_kf() const { return me_kf_; }
		void		   me_kf(const double &d) { me_kf_ = d; }
		const double & me_km() const { return me_km_; }
		void		   me_km(const double &d) { me_km_ = d; }

		// plastic potential parameters
		const double& me_mfg() const { return me_mfg_; }
		void		   me_mfg(const double& d) { me_mfg_ = d; }
		const double& me_alphag() const { return me_alphag_; }
		void		   me_alphag(const double& d) { me_alphag_ = d; }
		const double& me_betag() const { return me_betag_; }
		void		   me_betag(const double& d) { me_betag_ = d; }

		// wire sliding model
		const double& me_rata() const { return me_rata_; }
		void		   me_rata(const double& d) { me_rata_ = d; }
		const double& me_ratb() const { return me_ratb_; }
		void		   me_ratb(const double& d) { me_ratb_ = d; }

		// wire damage model
		const double& me_wpl1() const { return me_wpl1_; }
		void		   me_wpl1(const double& d) { me_wpl1_ = d; }
		const double& me_wpl2() const { return me_wpl2_; }
		void		   me_wpl2(const double& d) { me_wpl2_ = d; }
		const double& me_dalpha() const { return me_dalpha_; }
		void		   me_dalpha(const double& d) { me_dalpha_ = d; }

		// state variables
		const double & me_nf() const { return me_nf_; }
		void		   me_nf(const double &d) { me_nf_ = d; }
		const double & me_u_p() const { return me_u_p_; }
		void		   me_u_p(const double &d) { me_u_p_ = d; }
		const double& me_omega_p() const { return me_omega_p_; }
		void		   me_omega_p(const double& d) { me_omega_p_ = d; }
		const double & me_force() const { return me_force_; }
		void		   me_force(const double &d) { me_force_ = d; }
		const double & me_moment() const { return me_moment_; }
		void		   me_moment(const double &d) { me_moment_ = d; }
		const double& me_force_damaged() const { return me_force_damaged_; }
		void		   me_force_damaged(const double& d) { me_force_damaged_ = d; }
		const double& me_moment_damaged() const { return me_moment_damaged_; }
		void		   me_moment_damaged(const double& d) { me_moment_damaged_ = d; }
		const double& me_u_slide() const { return me_u_slide_; }
		void		   me_u_slide(const double& d) { me_u_slide_ = d; }
		const int& me_moment_sign() const { return me_moment_sign_; }
		void		   me_moment_sign(const int& d) { me_moment_sign_ = d; }
		const double& me_wp() const { return me_wp_; } // never gonna set this though
		void		   me_wp(const double& d) { me_wp_ = d; }
		const bool & me_isplastic() const { return me_isplastic_; }
		void		   me_isplastic(const bool &d) { me_isplastic_ = d; }

		bool     hasDamping() const { return dpProps_ ? true : false; }
		double   dp_nratio() const { return (hasDamping() ? (dpProps_->dp_nratio_) : 0.0); }
		void     dp_nratio(const double &d) { if (!hasDamping()) return; dpProps_->dp_nratio_ = d; }
		double   dp_sratio() const { return hasDamping() ? dpProps_->dp_sratio_ : 0.0; }
		void     dp_sratio(const double &d) { if (!hasDamping()) return; dpProps_->dp_sratio_ = d; }
		int      dp_mode() const { return hasDamping() ? dpProps_->dp_mode_ : -1; }
		void     dp_mode(int i) { if (!hasDamping()) return; dpProps_->dp_mode_ = i; }
		DVect    dp_F() const { return hasDamping() ? dpProps_->dp_F_ : DVect(0.0); }
		void     dp_F(const DVect &f) { if (!hasDamping()) return; dpProps_->dp_F_ = f; }

		bool    hasEnergies() const { return energies_ ? true : false; }
		double  estrain() const { return hasEnergies() ? energies_->estrain_ : 0.0; }
		void    estrain(const double &d) { if (!hasEnergies()) return; energies_->estrain_ = d; }
		double  eslip() const { return hasEnergies() ? energies_->eslip_ : 0.0; }
		void    eslip(const double &d) { if (!hasEnergies()) return; energies_->eslip_ = d; }
		double  edashpot() const { return hasEnergies() ? energies_->edashpot_ : 0.0; }
		void    edashpot(const double &d) { if (!hasEnergies()) return; energies_->edashpot_ = d; }
		double  epbstrain() const { return hasEnergies() ? energies_->epbstrain_ : 0.0; }
		void    epbstrain(const double &d) { if (!hasEnergies()) return; energies_->epbstrain_ = d; }


		double me_eelastic() const { return hasEnergies() ? energies_->me_estrain_ : 0.0; }
		void me_eelastic(const double &d) { if (!hasEnergies()) return; energies_->me_estrain_ = d; }
		double me_eplastic() const { return hasEnergies() ? energies_->me_eplastic_ : 0.0; }
		void me_eplastic(const double &d) { if (!hasEnergies()) return; energies_->me_eplastic_ = d; }
		double me_ebelastic() const { return hasEnergies() ? energies_->me_ebstrain_ : 0.0; }
		void me_ebelastic(const double &d) { if (!hasEnergies()) return; energies_->me_ebstrain_ = d; }
		double me_ebplastic() const { return hasEnergies() ? energies_->me_ebplastic_ : 0.0; }
		void me_ebplastic(const double &d) { if (!hasEnergies()) return; energies_->me_ebplastic_ = d; }


		bool     hasPBond() const { return pbProps_ ? true : false; }
		int      pbState() const { return hasPBond() ? pbProps_->pb_state_ : 0; }
		void     pbState(int i) { if (!hasPBond()) return; pbProps_->pb_state_ = i; }
		double   pbRmul() const { return (hasPBond() ? (pbProps_->pb_rmul_) : 0.0); }
		void     pbRmul(const double &d) { if (!hasPBond()) return; pbProps_->pb_rmul_ = d; }
		double   pbKn() const { return (hasPBond() ? (pbProps_->pb_kn_) : 0.0); }
		void     pbKn(const double &d) { if (!hasPBond()) return; pbProps_->pb_kn_ = d; }
		double   pbKs() const { return (hasPBond() ? (pbProps_->pb_ks_) : 0.0); }
		void     pbKs(const double &d) { if (!hasPBond()) return; pbProps_->pb_ks_ = d; }
		double   pbMCF() const { return (hasPBond() ? (pbProps_->pb_mcf_) : 0.0); }
		void     pbMCF(const double &d) { if (!hasPBond()) return; pbProps_->pb_mcf_ = d; }
		double   pbTen() const { return (hasPBond() ? (pbProps_->pb_ten_) : 0.0); }
		void     pbTen(const double &d) { if (!hasPBond()) return; pbProps_->pb_ten_ = d; }
		double   pbCoh() const { return (hasPBond() ? (pbProps_->pb_coh_) : 0.0); }
		void     pbCoh(const double &d) { if (!hasPBond()) return; pbProps_->pb_coh_ = d; }
		double   pbFA() const { return (hasPBond() ? (pbProps_->pb_fa_) : 0.0); }
		void     pbFA(const double &d) { if (!hasPBond()) return; pbProps_->pb_fa_ = d; }
		DVect    pbF() const { return hasPBond() ? pbProps_->pb_F_ : DVect(0.0); }
		void     pbF(const DVect &f) { if (!hasPBond()) return; pbProps_->pb_F_ = f; }
		DAVect   pbM() const { return hasPBond() ? pbProps_->pb_M_ : DAVect(0.0); }
		void     pbM(const DAVect &m) { if (!hasPBond()) return; pbProps_->pb_M_ = m; }
		DVect2   pbTransStiff() const { return hasPBond() ? pbProps_->pbTransStiff_ : DVect2(0.0); }
		void     pbTransStiff(const DVect2 &f) { if (!hasPBond()) return; pbProps_->pbTransStiff_ = f; }
		DAVect   pbAngStiff() const { return hasPBond() ? pbProps_->pbAngStiff_ : DAVect(0.0); }
		void     pbAngStiff(const DAVect &m) { if (!hasPBond()) return; pbProps_->pbAngStiff_ = m; }

		uint inheritanceField() const { return inheritanceField_; }
		void inheritanceField(uint i) { inheritanceField_ = i; }

		const DVect2 & effectiveTranslationalStiffness()  const { return effectiveTranslationalStiffness_; }
		void           effectiveTranslationalStiffness(const DVect2 &v) { effectiveTranslationalStiffness_ = v; }

		/// Return the total force that the contact model holds.
		virtual DVect    getForce(const IContactMechanical *) const;

		/// Return the total M on 1 that the contact model holds
		virtual DAVect   getMOn1(const IContactMechanical *) const;

		/// Return the total M on 1 that the contact model holds
		virtual DAVect   getMOn2(const IContactMechanical *) const;


		/// Return the total moment on 1 that the contact model holds
		virtual DAVect   getMomentOn1(const IContactMechanical*) const;

		/// Return the total moment on 1 that the contact model holds
		virtual DAVect   getMomentOn2(const IContactMechanical*) const;




	private:
		static int index_;

		struct Energies {
			Energies() : estrain_(0.0), eslip_(0.0), edashpot_(0.0), epbstrain_(0.0), me_estrain_(0.0), me_eplastic_(0.0), me_ebstrain_(0.0), me_ebplastic_(0.0) {}
			double estrain_;  // elastic energy stored in contact 
			double eslip_;    // work dissipated by friction 
			double edashpot_;    // work dissipated by dashpots
			double epbstrain_; // parallel bond strain energy
			double me_estrain_;
			double me_eplastic_;
			double me_ebstrain_;
			double me_ebplastic_;
		};

		struct meEnergies {
			meEnergies() : me_estrain_(0.0),me_eplastic_(0.0),me_ebstrain_(0.0),me_ebplastic_(0.0){}
			double me_estrain_;
			double me_eplastic_;
			double me_ebstrain_;
			double me_ebplastic_;
		};

		struct dpProps {
			dpProps() : dp_nratio_(0.0), dp_sratio_(0.0), dp_mode_(0), dp_F_(DVect(0.0)) {}
			double dp_nratio_;    // normal viscous critical damping ratio
			double dp_sratio_;    // shear  viscous critical damping ratio
			int    dp_mode_;      // for viscous mode (0-4) 0 = dashpots, 1 = tensile limit, 2 = shear limit, 3 = limit both
			DVect  dp_F_;		  // Force in the dashpots
		};

		struct pbProps {
			pbProps() : pb_state_(0), pb_rmul_(1.0), pb_kn_(0.0), pb_ks_(0.0),
				pb_mcf_(1.0), pb_ten_(0.0), pb_coh_(0.0), pb_fa_(0.0), pb_F_(DVect(0.0)), pb_M_(DAVect(0.0)),
				pbTransStiff_(0.0), pbAngStiff_(0.0) {}
			// parallel bond
			int     pb_state_;        // Bond mode - 0 (NBNF), 1 (NBFT), 2 (NBFS), 3 (B)
			double  pb_rmul_;         // Radius multiplier
			double  pb_kn_;           // normal stiffness
			double  pb_ks_;           // shear stiffness
			double  pb_mcf_;          // M contribution factor 
			double  pb_ten_;          // normal strength 
			double  pb_coh_;          // cohesion
			double  pb_fa_;           // friction angle
			DVect   pb_F_;            // Force in parallel bond
			DAVect  pb_M_;            // M in parallel bond
			DVect2  pbTransStiff_;    // (Normal,Shear) Translational stiffness of the parallel bond
			DAVect  pbAngStiff_;      // (Normal,Shear) Rotational    stiffness of the parallel bond
		};

		bool   updateKn(const IContactMechanical *con);
		bool   updateKs(const IContactMechanical *con);
		bool   updateFric(const IContactMechanical *con);
		double pbStrainEnergy() const;                     // Compute  bond strain energy 

		void   updateEffectiveStiffness(ContactModelMechanicalState *state);

		DVect3 pbData(const IContactMechanical *con) const; // Bond area and inertia
		DVect2 pbSMax(const IContactMechanical *con) const; // Maximum stress (tensile,shear) at bond periphery
		double pbShearStrength(const double &pbArea) const; // Bond shear strength
		void   setDampCoefficients(const double &mass, double *vcn, double *vcs);

		// inheritance fields
		quint32 inheritanceField_;

		// linear model
		double      kn_;        // normal stiffness
		double      ks_;        // shear stiffness
		double      fric_;      // Coulomb friction coefficient
		DVect       lin_F_;     // Force carried in the linear model
		bool        lin_S_;     // the current sliding state
		uint        lin_mode_;  // specifies absolute (0) or incremental (1) behavior for the the linear part 
		double      rgap_;      // reference gap for the linear part

		dpProps *   dpProps_;     // The viscous properties
		pbProps *   pbProps_;     // The parallel bond properties

		double      userArea_;    // User specified area 

		Energies *   energies_;   // energies
		DVect2  effectiveTranslationalStiffness_;



		double me_ea_;	     // axial stiffness
		double me_ei_;		 // bending stiffness
		double me_d_;	     // wire diameter. this is not part of the macro-element model itself but is used to link it with shear/twisting  
		double me_ref_l_;	 // reference wire length
		double me_curr_l_;   // current wire length

		// yield surface
		double me_nzero_;		 // material's yield stress 
		double me_nult_;		 // material's maximum stress (never reached)
		double me_ntens_;		 // buckling limit (empirical)

		double me_mff_;			 // shape parameters
		double me_betaf_;
		
		// plastic potential
		double me_mfg_;			 // shape parameters
		double me_alphag_;
		double me_betag_;

		// hardening
		double me_kf_;		 // hardening parameter for pure tensile plastic deformation
		double me_km_;		 // hardening parameter for pure bending plastic deformation
		
		// damage
		double me_wpl1_;	// damage onset
		double me_wpl2_;	// damage end (failure)
		double me_dalpha_;	// damage rate

		// sliding
		double me_rata_;	// sliding rate (friction + angle)
		double me_ratb_;	// maximum sliding displacement

		// state variables
		double me_nf_;
		double me_u_p_;
		double me_omega_p_;
		double me_force_;
		double me_force_damaged_;
		double me_moment_;
		double me_moment_damaged_;
		int	   me_moment_sign_;
		double me_u_slide_; //for funzies
		double me_wp_;
		bool   me_isplastic_;
		

	};
} // namespace itascaxd
// EoF