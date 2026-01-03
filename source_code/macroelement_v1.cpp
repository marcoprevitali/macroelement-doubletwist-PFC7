/** \file contactmodelexample.cpp
  * \addtogroup ContactModelExample
  * @{
*/
#include "macroelement_v1.h"
#include "contactmodel/src/contactmodelthermal.h"
#include "fish/src/parameter.h"
#include "utility/src/tptr.h"
#include "shared/src/mathutil.h"
#include "version.txt"

#include "kernel/interface/iprogram.h"
#include "module/interface/icontact.h"
#include "module/interface/icontactmechanical.h"
#include "module/interface/icontactthermal.h"
#include "module/interface/ifishcalllist.h"
#include "module/interface/ipiecemechanical.h"
#include "module/interface/ipiece.h"


#ifdef macroelement_LIB
int __stdcall DllMain(void *, unsigned, void *) {
	return 1;
}

extern "C" EXPORT_TAG const char *getName() {
#if DIM==3
	return "contactmodelmechanical3dmacroelement";
#else
	return "contactmodelmechanical2dmacroelement";
#endif
}

extern "C" EXPORT_TAG unsigned getMajorVersion() {
	return MAJOR_VERSION;
}
extern "C" EXPORT_TAG unsigned getMinorVersion() {
	return UPDATE_VERSION;
}


extern "C" EXPORT_TAG void *createInstance() {
	cmodelsxd::ContactModelmacroelement *m = new cmodelsxd::ContactModelmacroelement();
	return (void *)m;
}
#endif 
namespace cmodelsxd {
	static const quint32 linKnMask = 0x00002; // Base 1!
	static const quint32 linKsMask = 0x00004;
	static const quint32 linFricMask = 0x00008;

	using namespace itasca;

	int ContactModelmacroelement::index_ = -1;
	UInt ContactModelmacroelement::getMinorVersion() const { return MINOR_VERSION; }

	ContactModelmacroelement::ContactModelmacroelement() : inheritanceField_(linKnMask | linKsMask | linFricMask)
		, kn_(0.0)
		, ks_(0.0) 
		, fric_(0.0)
		, me_ea_(0.0) 
		, me_ei_(0.0)
		, me_d_(0.0027)
		, me_ref_l_(0.06)
		, me_curr_l_(me_ref_l_)
		, me_nzero_(0.0)
		, me_nult_(0.0)
		, me_ntens_(0.0) 
		, me_mff_(0.0)
		, me_betaf_(0.0)
		, me_kf_(0.0)
		, me_km_(0.0)
		, me_mfg_(0.0)
		, me_alphag_(0.0)
		, me_betag_(0.0)
		, me_wpl1_(0.0)
		, me_wpl2_(0.0)
		, me_dalpha_(0.0)
		, me_rata_(0.0)
		, me_ratb_(0.0)
		, me_nf_(me_nzero_)
		, me_u_p_(0.0)
		, me_omega_p_(0.0)
		, me_force_(0.0)
		, me_moment_(0.0)
		, me_isplastic_(false)
		, me_u_slide_(0.0)
		, me_wp_(0.0)
		, me_force_damaged_(0.0)
		, me_moment_damaged_(0.0)
		, me_moment_sign_(1)
		, lin_F_(DVect(0.0))
		, lin_S_(false)
		, lin_mode_(0)
		, rgap_(0.0)
		, dpProps_(0)
		, pbProps_(0)
		, userArea_(0)
		, energies_(0)
		, effectiveTranslationalStiffness_(DVect2(0.0)) {
		//    setFromParent(ContactModelMechanicalList::instance()->find(getName()));
	}

	ContactModelmacroelement::~ContactModelmacroelement() {
		if (dpProps_)
			delete dpProps_;
		if (pbProps_)
			delete pbProps_;
		if (energies_)
			delete energies_;
	}

	void ContactModelmacroelement::archive(ArchiveStream &stream) {
		stream & kn_;
		stream & ks_;
		stream & fric_;
		stream & me_ea_;
		stream & me_ei_;
		stream& me_d_;
		stream& me_ref_l_;
		stream& me_curr_l_;
		stream& me_nzero_;
		stream& me_nult_;
		stream& me_ntens_;
		stream& me_mff_;
		stream& me_betaf_;
		stream& me_kf_;
		stream& me_km_;
		stream& me_mfg_;
		stream& me_alphag_;
		stream& me_betag_;
		stream& me_wpl1_;
		stream& me_wpl2_;
		stream& me_dalpha_;
		stream& me_rata_;
		stream& me_ratb_;
		stream & me_nf_;
		stream& me_u_p_;
		stream & me_omega_p_;
		stream & me_force_;
		stream & me_moment_;
		stream& me_moment_sign_;
		stream& me_u_slide_;
		stream& me_isplastic_;
		stream& me_wp_;
		stream& me_force_damaged_;
		stream& me_moment_damaged_;
		stream & lin_F_;
		stream & lin_S_;
		stream & lin_mode_;
		if (stream.getArchiveState() == ArchiveStream::Save) {
			bool b = false;
			if (dpProps_) {
				b = true;
				stream & b;
				stream & dpProps_->dp_nratio_;
				stream & dpProps_->dp_sratio_;
				stream & dpProps_->dp_mode_;
				stream & dpProps_->dp_F_;
			}
			else
				stream & b;

			b = false;
			if (energies_) {
				b = true;
				stream & b;
				stream & energies_->estrain_;
				stream & energies_->eslip_;
				stream & energies_->edashpot_;
				stream & energies_->epbstrain_;
				
			}
			else
				stream & b;
			b = false;
			if (pbProps_) {
				b = true;
				stream & b;
				stream & pbProps_->pb_state_;
				stream & pbProps_->pb_rmul_;
				stream & pbProps_->pb_kn_;
				stream & pbProps_->pb_ks_;
				stream & pbProps_->pb_mcf_;
				stream & pbProps_->pb_ten_;
				stream & pbProps_->pb_coh_;
				stream & pbProps_->pb_fa_;
				stream & pbProps_->pb_F_;
				stream & pbProps_->pb_M_;
			}
			else
				stream & b;
		}
		else {
			bool b(false);
			stream & b;
			if (b) {
				if (!dpProps_)
					dpProps_ = NEWC(dpProps());
				stream & dpProps_->dp_nratio_;
				stream & dpProps_->dp_sratio_;
				stream & dpProps_->dp_mode_;
				stream & dpProps_->dp_F_;
			}
			stream & b;
			if (b) {
				if (!energies_)
					energies_ = NEWC(Energies());
				stream & energies_->estrain_;
				stream & energies_->eslip_;
				stream & energies_->edashpot_;
				stream & energies_->epbstrain_;
				stream & energies_->me_estrain_;
				stream & energies_->me_eplastic_;
				stream & energies_->me_ebstrain_;
				stream & energies_->me_ebplastic_;
			}

			stream & b;
			if (b) {
				if (!pbProps_)
					pbProps_ = NEWC(pbProps());
				stream & pbProps_->pb_state_;
				stream & pbProps_->pb_rmul_;
				stream & pbProps_->pb_kn_;
				stream & pbProps_->pb_ks_;
				stream & pbProps_->pb_mcf_;
				stream & pbProps_->pb_ten_;
				stream & pbProps_->pb_coh_;
				stream & pbProps_->pb_fa_;
				stream & pbProps_->pb_F_;
				stream & pbProps_->pb_M_;
			}
		}

		stream & inheritanceField_;
		stream & effectiveTranslationalStiffness_;

		if (stream.getArchiveState() == ArchiveStream::Save || stream.getRestoreVersion() == getMinorVersion())
			stream & rgap_;

		if (stream.getArchiveState() == ArchiveStream::Save || stream.getRestoreVersion() > 1)
			stream & userArea_;
	}

	void ContactModelmacroelement::copy(const ContactModel *cm) {
		ContactModelMechanical::copy(cm);
		const ContactModelmacroelement *in = dynamic_cast<const ContactModelmacroelement*>(cm);
		if (!in) throw std::runtime_error("Internal error: contact model dynamic cast failed.");
		kn(in->kn());
		ks(in->ks());
		fric(in->fric());
		me_ea(in->me_ea());
		me_ei(in->me_ei());
		me_d(in->me_d());
		me_ref_l(in->me_ref_l());
		me_curr_l(in->me_curr_l());
		me_nzero(in->me_nzero());
		me_nult(in->me_nult());
		me_ntens(in->me_ntens());
		me_mff(in->me_mff());
		me_betaf(in->me_betaf());
		me_mfg(in->me_mfg());
		me_alphag(in->me_alphag());
		me_betag(in->me_betag());
		me_kf(in->me_kf());
		me_km(in->me_km());
		me_wpl1(in->me_wpl1());
		me_wpl2(in->me_wpl2());
		me_dalpha(in->me_dalpha());
		me_rata(in->me_rata());
		me_ratb(in->me_ratb());
		me_nf(in->me_nf());
		me_u_p(in->me_u_p());
		me_omega_p(in->me_omega_p());
		me_force(in->me_force());
		me_moment(in->me_moment());
		me_isplastic(in->me_isplastic());
		me_moment_sign(in->me_moment_sign());
		me_u_slide(in->me_u_slide());
		me_wp(in->me_wp());
		me_force_damaged(in->me_force_damaged());
		me_moment_damaged(in->me_moment_damaged());
		lin_F(in->lin_F());
		lin_S(in->lin_S());
		lin_mode(in->lin_mode());
		rgap(in->rgap());
		if (in->hasDamping()) {
			if (!dpProps_)
				dpProps_ = NEWC(dpProps());
			dp_nratio(in->dp_nratio());
			dp_sratio(in->dp_sratio());
			dp_mode(in->dp_mode());
			dp_F(in->dp_F());
		}
		if (in->hasEnergies()) {
			if (!energies_)
				energies_ = NEWC(Energies());
			estrain(in->estrain());
			eslip(in->eslip());
			edashpot(in->edashpot());
			epbstrain(in->epbstrain());
			me_eelastic(in->me_eelastic());
			me_eplastic(in->me_eplastic());
			me_ebelastic(in->me_ebelastic());
			me_ebplastic(in->me_ebplastic());
		}
		if (in->hasPBond()) {
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbState(in->pbState());
			pbRmul(in->pbRmul());
			pbKn(in->pbKn());
			pbKs(in->pbKs());
			pbMCF(in->pbMCF());
			pbTen(in->pbTen());
			pbCoh(in->pbCoh());
			pbFA(in->pbFA());
			pbF(in->pbF());
			pbM(in->pbM());
			pbTransStiff(in->pbTransStiff());
			pbAngStiff(in->pbAngStiff());
		}
		userArea_ = in->userArea_;
		inheritanceField(in->inheritanceField());
		effectiveTranslationalStiffness(in->effectiveTranslationalStiffness());
	}

	QVariant ContactModelmacroelement::getProperty(uint i, const IContact *con) const {
		QVariant var;
		switch (i) {
		case kwLinKn:        return kn_;
		case kwLinKs:        return ks_;
		case kwLinFric:      return fric_;
		case kwEA:			 return me_ea_;
		case kwEI:			 return me_ei_;
		case kwD:			 return me_d_;
		case kwRefL:		 return me_ref_l_;
		case kwCurrL:		 return me_curr_l_;
		case kwNZero:		 return me_nzero_;
		case kwNult:		 return me_nult_;
		case kwNtens:		 return me_ntens_;
		case kwMff:			 return me_mff_;
		case kwBetaf:		 return me_betaf_;
		case kwKF:			 return me_kf_;
		case kwKM:			 return me_km_;
		case kwMfg:			 return me_mfg_;
		case kwAlphag:		 return me_alphag_;
		case kwBetag:		 return me_betag_;
		case kwWpl1:		 return me_wpl1_;
		case kwWpl2:		 return me_wpl2_;
		case kwDalpha:		 return me_dalpha_;
		case kwRatA:		 return me_rata_;
		case kwRatB:		 return me_ratb_;
		case kwNf:			 return me_nf_;
		case kwUP:			 return me_u_p_;
		case kwOmegaP:		 return me_omega_p_;
		case kwForce:		 return me_force_;
		case kwMoment:		 return me_moment_;
		case kwIsPlastic:	 return me_isplastic_;
		case kwUSlide:		 return me_u_slide_;
		case kwWP:			 return me_wp_;
		case kwForceDamaged: return me_force_damaged_;
		case kwMomentDamaged: return me_moment_damaged_;
		case kwMomentSign:	return me_moment_sign_;
		case kwLinMode:      return lin_mode_;
		case kwLinF:         var.setValue(lin_F_); return var;
		case kwLinS:		 return lin_S_;
		case kwRGap:	     return rgap_;
		case kwEmod: {
			const IContactMechanical *c(convert_getcast<IContactMechanical>(con));
			if (c == nullptr) return 0.0;
			double rsq(std::max(c->getEnd1Curvature().y(), c->getEnd2Curvature().y()));
			double rsum(0.0);
			if (c->getEnd1Curvature().y())
				rsum += 1.0 / c->getEnd1Curvature().y();
			if (c->getEnd2Curvature().y())
				rsum += 1.0 / c->getEnd2Curvature().y();
			if (userArea_) {
#ifdef THREED
				rsq = std::sqrt(userArea_ / dPi);
#else
				rsq = userArea_ / 2.0;
#endif        
				rsum = rsq + rsq;
				rsq = 1. / rsq;
			}
#ifdef TWOD             
			return (kn_ * rsum * rsq / 2.0);
#else                     
			return (kn_ * rsum * rsq * rsq) / dPi;
#endif                    
		}
		case kwKRatio:      return (ks_ == 0.0) ? 0.0 : (kn_ / ks_);
		case kwDpNRatio:    return dpProps_ ? dpProps_->dp_nratio_ : 0;
		case kwDpSRatio:    return dpProps_ ? dpProps_->dp_sratio_ : 0;
		case kwDpMode:      return dpProps_ ? dpProps_->dp_mode_ : 0;
		case kwUserArea:    return userArea_;
		case kwDpF: {
			dpProps_ ? var.setValue(dpProps_->dp_F_) : var.setValue(DVect(0.0));
			return var;
		}
		case kwPbState:     return pbProps_ ? pbProps_->pb_state_ : 0;
		case kwPbRMul:      return pbProps_ ? pbProps_->pb_rmul_ : 1.0;
		case kwPbKn:        return pbProps_ ? pbProps_->pb_kn_ : 0;
		case kwPbKs:        return pbProps_ ? pbProps_->pb_ks_ : 0;
		case kwPbMcf:       return pbProps_ ? pbProps_->pb_mcf_ : 1.0;
		case kwPbTStrength: return pbProps_ ? pbProps_->pb_ten_ : 0.0;
		case kwPbSStrength: {
			if (!pbProps_) return 0.0;
			const IContactMechanical *c(convert_getcast<IContactMechanical>(con));
			double pbArea = pbData(c).x();
			return pbShearStrength(pbArea);
		}
		case kwPbCoh:       return pbProps_ ? pbProps_->pb_coh_ : 0;
		case kwPbFa:        return pbProps_ ? pbProps_->pb_fa_ : 0;
		case kwPbSig: {
			if (!pbProps_ || pbProps_->pb_state_ < 3) return 0.0;
			const IContactMechanical *c(convert_getcast<IContactMechanical>(con));
			return pbSMax(c).x();
		}
		case kwPbTau: {
			if (!pbProps_ || pbProps_->pb_state_ < 3) return 0.0;
			const IContactMechanical *c(convert_getcast<IContactMechanical>(con));
			return pbSMax(c).y();
		}
		case kwPbF: {
			pbProps_ ? var.setValue(pbProps_->pb_F_) : var.setValue(DVect(0.0));
			return var;
		}
		case kwPbM: {
			pbProps_ ? var.setValue(pbProps_->pb_M_) : var.setValue(DAVect(0.0));
			return var;
		}
		case kwPbRadius: {
			if (!pbProps_) return 0.0;
			const IContactMechanical *c(convert_getcast<IContactMechanical>(con));
			double Cmax1 = c->getEnd1Curvature().y();
			double Cmax2 = c->getEnd2Curvature().y();
			double br = pbProps_->pb_rmul_ * 1.0 / std::max(Cmax1, Cmax2);
			if (userArea_)
#ifdef THREED
				br = std::sqrt(userArea_ / dPi);
#else
				br = userArea_ / 2.0;
#endif
			return br;
		}
		case kwPbEmod: {
			if (!pbProps_) return 0.0;
			const IContactMechanical *c(convert_getcast<IContactMechanical>(con));
			double rsum(0.0);
			if (c->getEnd1Curvature().y())
				rsum += 1.0 / c->getEnd1Curvature().y();
			if (c->getEnd2Curvature().y())
				rsum += 1.0 / c->getEnd2Curvature().y();
			if (userArea_) {
#ifdef THREED
				double rad = std::sqrt(userArea_ / dPi);
#else
				double rad = userArea_ / 2.0;
#endif        
				rsum = 2 * rad;
			}
			double emod = pbProps_->pb_kn_ * rsum;
			return emod;
		}
		case kwPbKRatio: {
			if (!pbProps_) return 0.0;
			return (pbProps_->pb_ks_ == 0.0) ? 0.0 : (pbProps_->pb_kn_ / pbProps_->pb_ks_);
		}
		}
		//assert(0); // technically this should be a way to check if there are empty keywords or something so this might cause bugs if left commented
		return QVariant();
	}

	bool ContactModelmacroelement::getPropertyGlobal(uint i) const {
		switch (i) {
		case kwLinF:
		case kwDpF:
		case kwPbF:
			return false;
		}
		return true;
	}

	bool ContactModelmacroelement::setProperty(uint i, const QVariant &v, IContact *) {
		pbProps pb;
		dpProps dp;

		switch (i) {
		case kwLinKn: {
			if (!v.canConvert<double>())
				throw Exception("kn must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative kn not allowed.");
			kn_ = val;
			return true;
		}
		case kwLinKs: {
			if (!v.canConvert<double>())
				throw Exception("ks must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative ks not allowed.");
			ks_ = val;
			return true;
		}
		case kwLinFric: {
			if (!v.canConvert<double>())
				throw Exception("fric must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative fric not allowed.");
			fric_ = val;
			return false;
		}

		case kwEA: {
			if (!v.canConvert<double>())
				throw Exception("Axial Stiffness EA must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Axial Stiffness EA must be a positive value.");
			me_ea_ = val;
			kn_ = me_ea_ / 0.06; // for automatic step calc
			double me_poisson_ = 0.3;
			ks_ = kn_ * 0.5 / (1.0 + me_poisson_);
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->pb_kn_ = kn_ / ((me_d_ / 2) * (me_d_ / 2) * M_PI);
			pbProps_->pb_ks_ = ks_ / ((me_d_ / 2) * (me_d_ / 2) * M_PI);
			return false;
		}
		case kwEI: {
			if (!v.canConvert<double>())
				throw Exception("Flexural stiffness EI must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Flexural stiffness EI must be a positive value.");
			me_ei_ = val;
			return false;
		}
		case kwD: {
			if (!v.canConvert<double>())
				throw Exception("Wire diameter D must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Wire diameter D must be a positive value.");
			me_d_ = val;
			return false;
		}
		case kwRefL: {
			if (!v.canConvert<double>())
				throw Exception("Reference Wire length L must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Reference Wire length L must be a positive value.");
			me_ref_l_ = val;
			return false;
		}
		case kwCurrL: {
			if (!v.canConvert<double>())
				throw Exception("Current Wire length L must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Current Wire length L must be a positive value.");
			me_curr_l_ = val;
			return false;
		}
		case kwNZero: {
			if (!v.canConvert<double>())
				throw Exception("Yield stress must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Yield stress must be a positive value.");
			me_nzero_ = val;
			return false;
		}
		case kwNult: {
			if (!v.canConvert<double>())
				throw Exception("Young modulus must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Young modulus must be a positive value.");
			me_nult_ = val;
			return false;
		}
		case kwNtens: {
			if (!v.canConvert<double>())
				throw Exception("Buckling limit must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Buckling limit must be a positive value.");
			me_ntens_ = val;
			return false;
		}
		case kwMff: {
			if (!v.canConvert<double>())
				throw Exception("Yield function Mf must be a double");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Yield function Mf must be a positive value.");
			me_mff_ = val;
			return false;
		}
		case kwBetaf: {
			if (!v.canConvert<double>())
				throw Exception("Yield function Beta parameter must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Yield function Beta must be a positive value.");
			me_betaf_ = val;
			return false;
		}
		case kwKF: {
			if (!v.canConvert<double>())
				throw Exception("KF must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("KF must be a positive value.");
			me_kf_ = val;
			return false;
		}
		case kwKM: {
			if (!v.canConvert<double>())
				throw Exception("KM must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("KM must be a positive value.");
			me_km_ = val;
			return false;
		}
		case kwMfg: {
			if (!v.canConvert<double>())
				throw Exception("Plastic potential Mf parameter must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Plastic potential Mf parameter must be a positive value.");
			me_mfg_ = val;
			return false;
		}
		case kwAlphag: {
			if (!v.canConvert<double>())
				throw Exception("Plastic potential Alpha parameter must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Plastic potential Alpha parameter must be a positive value.");
			me_alphag_ = val;
			return false;
		}
		case kwBetag: {
			if (!v.canConvert<double>())
				throw Exception("Plastic potential Beta parameter must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Plastic potential Beta parameter must be a positive value.");
			me_betag_ = val;
			return false;
		}
		case kwWpl1: {
			if (!v.canConvert<double>())
				throw Exception("Wire damage Wpl1 parameter must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Wire damage Wpl1 parameter must be a positive value.");
			me_wpl1_ = val;
			return false;
					}
		case kwWpl2: {
			if (!v.canConvert<double>())
				throw Exception("Wire damage Wpl2 parameter must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Wire damage Wpl2 parameter must be a positive value.");
			me_wpl2_ = val;
			return false;
		}
		case kwDalpha: {
			if (!v.canConvert<double>())
				throw Exception("Wire damage D alpha parameter must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Wire damage D alpha parameter must be a positive value.");
			me_dalpha_ = val;
			return false;
		}
		case kwRatA: {
			if (!v.canConvert<double>())
				throw Exception("Nonlinear wire-sliding elasticity a parameter must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Nonlinear wire-sliding elasticity a parameter must be a positive value.");
			me_rata_ = val;
			return false;
		}
		case kwRatB: {
			if (!v.canConvert<double>())
				throw Exception("Nonlinear wire-sliding elasticity u_distortion parameter must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Nonlinear wire-sliding elasticity u_distortion parameter must be a positive value.");
			me_ratb_ = val;
			return false;
		}
		case kwLinF: {
			if (!v.canConvert<DVect>())
				throw Exception("lin_force must be a vector.");
			DVect val(v.value<DVect>());
			lin_F_ = val;
			return false;
		}
		case kwLinMode: {
			if (!v.canConvert<uint>())
				throw Exception("lin_mode must be 0 (absolute) or 1 (incremental).");
			uint val(v.toUInt());
			if (val > 1)
				throw Exception("lin_mode must be 0 (absolute) or 1 (incremental).");
			lin_mode_ = val;
			return false;
		}
		case kwRGap: {
			if (!v.canConvert<double>())
				throw Exception("Reference gap must be a double.");
			double val(v.toDouble());
			rgap_ = val;
			return false;
		}
		case kwPbRMul: {
			if (!v.canConvert<double>())
				throw Exception("pb_rmul must be a double.");
			double val(v.toDouble());
			if (val <= 0.0)
				throw Exception("pb_rmul must be positive.");
			if (val == 1.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->pb_rmul_ = val;
			return false;
		}
		case kwPbKn: {
			if (!v.canConvert<double>())
				throw Exception("pb_kn must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative pb_kn not allowed.");
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->pb_kn_ = val;
			return true;
		}
		case kwPbKs: {
			if (!v.canConvert<double>())
				throw Exception("pb_ks must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative pb_ks not allowed.");
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->pb_ks_ = val;
			return true;
		}
		case kwPbMcf: {
			if (!v.canConvert<double>())
				throw Exception("pb_mcf must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative pb_mcf not allowed.");
			if (val > 1.0)
				throw Exception("pb_mcf must be lower or equal to 1.0.");
			if (val == 1.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->pb_mcf_ = val;
			return false;
		}
		case kwPbTStrength: {
			if (!v.canConvert<double>())
				throw Exception("pb_ten must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative pb_ten not allowed.");
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->pb_ten_ = val;
			return false;
		}
		case kwPbCoh: {
			if (!v.canConvert<double>())
				throw Exception("pb_coh must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative pb_coh not allowed.");
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->pb_coh_ = val;
			return false;
		}
		case kwPbFa: {
			if (!v.canConvert<double>())
				throw Exception("pb_fa must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative pb_fa not allowed.");
			if (val >= 90.0)
				throw Exception("pb_fa must be lower than 90.0 degrees.");
			if (val == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->pb_fa_ = val;
			return false;
		}
		case kwPbF: {
			if (!v.canConvert<DVect>())
				throw Exception("pb_force must be a vector.");
			DVect val(v.value<DVect>());
			if (val.fsum() == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->pb_F_ = val;
			return false;
		}
		case kwPbM: {
			DAVect val(0.0);
#ifdef TWOD               
			if (!v.canConvert<DAVect>() && !v.canConvert<double>())
				throw Exception("pb_M must be an angular vector.");
			if (v.canConvert<DAVect>())
				val = DAVect(v.value<DAVect>());
			else
				val = DAVect(v.toDouble());
#else
			if (!v.canConvert<DAVect>() && !v.canConvert<DVect>())
				throw Exception("pb_M must be an angular vector.");
			if (v.canConvert<DAVect>())
				val = DAVect(v.value<DAVect>());
			else
				val = DAVect(v.value<DVect>());
#endif
			if (val.fsum() == 0.0 && !pbProps_)
				return false;
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			pbProps_->pb_M_ = val;
			return false;
		}
		case kwDpNRatio: {
			if (!v.canConvert<double>())
				throw Exception("dp_nratio must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative dp_nratio not allowed.");
			if (val == 0.0 && !dpProps_)
				return false;
			if (!dpProps_)
				dpProps_ = NEWC(dpProps());
			dpProps_->dp_nratio_ = val;
			return true;
		}
		case kwDpSRatio: {
			if (!v.canConvert<double>())
				throw Exception("dp_sratio must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative dp_sratio not allowed.");
			if (val == 0.0 && !dpProps_)
				return false;
			if (!dpProps_)
				dpProps_ = NEWC(dpProps());
			dpProps_->dp_sratio_ = val;
			return true;
		}
		case kwDpMode: {
			if (!v.canConvert<int>())
				throw Exception("The viscous mode dp_mode must be 0, 1, 2, or 3.");
			int val(v.toInt());
			if (val == 0 && !dpProps_)
				return false;
			if (val < 0 || val > 3)
				throw Exception("The viscous mode dp_mode must be 0, 1, 2, or 3.");
			if (!dpProps_)
				dpProps_ = NEWC(dpProps());
			dpProps_->dp_mode_ = val;
			return false;
		}
		case kwDpF: {
			if (!v.canConvert<DVect>())
				throw Exception("dp_force must be a vector.");
			DVect val(v.value<DVect>());
			if (val.fsum() == 0.0 && !dpProps_)
				return false;
			if (!dpProps_)
				dpProps_ = NEWC(dpProps());
			dpProps_->dp_F_ = val;
			return false;
		}
		case kwUserArea: {
			if (!v.canConvert<double>())
				throw Exception("user_area must be a double.");
			double val(v.toDouble());
			if (val < 0.0)
				throw Exception("Negative user_area not allowed.");
			userArea_ = val;
			return true;
		}
		}
		//    assert(0);
		return false;
	}

	bool ContactModelmacroelement::getPropertyReadOnly(uint i) const {
		switch (i) {
		case kwDpF:
		case kwLinS:
		case kwEmod:
		case kwKRatio:
		case kwPbState:
		case kwPbRadius:
		case kwPbSStrength:
		case kwPbSig:
		case kwPbTau:
		case kwPbEmod:
		case kwPbKRatio:
		case kwNf:
		case kwUP:
		case kwCurrL:
		case kwOmegaP:
		case kwForce:
		case kwIsPlastic:
		case kwUSlide:
		case kwWP:
		case kwMomentSign:
		case kwForceDamaged:
		case kwMomentDamaged:
			return true;
		default:
			break;
		}
		return false;
	}

	bool ContactModelmacroelement::supportsInheritance(uint i) const {
		switch (i) {
		case kwLinKn:
		case kwLinKs:
		case kwEA:
		case kwEI:
		case kwLinFric:
		case kwD:
		case kwRefL:
		case kwNZero:
		case kwNult:
		case kwNtens:
		case kwMff:
		case kwBetaf:
		case kwKF:
		case kwKM:
		case kwMfg:
		case kwAlphag:
		case kwWpl1:
		case kwWpl2:
		case kwDalpha:
		case kwRatA:
		case kwRatB:
			return true;
		default:
			break;
		}
		return false;
	}

	QString  ContactModelmacroelement::getMethodArguments(uint i) const {
		QString def = QString();
		switch (i) {
		case kwDeformability:
			return "emod,kratio";
		case kwPbDeformability:
			return "emod,kratio";
		case kwPbBond:
			return "gap";
		case kwPbUnbond:
			return "gap";
		}
		return def;
	}

	bool ContactModelmacroelement::setMethod(uint i, const QVector<QVariant> &vl, IContact *con) {
		IContactMechanical *c(convert_getcast<IContactMechanical>(con));
		switch (i) {
		case kwDeformability: {
			double emod;
			double krat;
			if (vl.at(0).isNull())
				throw Exception("Argument emod must be specified with method deformability in contact model %1.", getName());
			emod = vl.at(0).toDouble();
			if (emod < 0.0)
				throw Exception("Negative emod not allowed in contact model %1.", getName());
			if (vl.at(1).isNull())
				throw Exception("Argument kratio must be specified with method deformability in contact model %1.", getName());
			krat = vl.at(1).toDouble();
			if (krat < 0.0)
				throw Exception("Negative linear stiffness ratio not allowed in contact model %1.", getName());
			double rsq(std::max(c->getEnd1Curvature().y(), c->getEnd2Curvature().y()));
			double rsum(0.0);
			if (c->getEnd1Curvature().y())
				rsum += 1.0 / c->getEnd1Curvature().y();
			if (c->getEnd2Curvature().y())
				rsum += 1.0 / c->getEnd2Curvature().y();
			if (userArea_) {
#ifdef THREED
				rsq = std::sqrt(userArea_ / dPi);
#else
				rsq = userArea_ / 2.0;
#endif        
				rsum = rsq + rsq;
				rsq = 1. / rsq;
			}
#ifdef TWOD
			kn_ = 2.0 * emod / (rsq * rsum);
#else
			kn_ = dPi * emod / (rsq * rsq * rsum);
#endif
			ks_ = (krat == 0.0) ? 0.0 : kn_ / krat;
			setInheritance(1, false);
			setInheritance(2, false);
			return true;
		}
		case kwPbDeformability: {
			//if (!pbProps_ || pbProps_->pb_state_ != 3) return false;
			double emod;
			double krat;
			if (vl.at(0).isNull())
				throw Exception("Argument emod must be specified with method pb_deformability in contact model %1.", getName());
			emod = vl.at(0).toDouble();
			if (emod < 0.0)
				throw Exception("Negative emod not allowed in contact model %1.", getName());
			if (vl.at(1).isNull())
				throw Exception("Argument kratio must be specified with method pb_deformability in contact model %1.", getName());
			krat = vl.at(1).toDouble();
			if (krat < 0.0)
				throw Exception("Negative parallel bond stiffness ratio not allowed in contact model %1.", getName());
			double rsum(0.0);
			if (c->getEnd1Curvature().y())
				rsum += 1.0 / c->getEnd1Curvature().y();
			if (c->getEnd2Curvature().y())
				rsum += 1.0 / c->getEnd2Curvature().y();
			if (!pbProps_)
				pbProps_ = NEWC(pbProps());
			if (userArea_)
#ifdef THREED
				rsum = 2 * std::sqrt(userArea_ / dPi);
#else
				rsum = 2 * userArea_ / 2.0;
#endif
			pbProps_->pb_kn_ = emod / rsum;
			pbProps_->pb_ks_ = (krat == 0.0) ? 0.0 : pbProps_->pb_kn_ / krat;
			return true;
		}
		case kwPbBond: {
			if (pbProps_ && pbProps_->pb_state_ == 3) return false;
			double mingap = -1.0 * limits<double>::max();
			double maxgap = 0;
			if (vl.at(0).canConvert<Double>())
				maxgap = vl.at(0).toDouble();
			else if (vl.at(0).canConvert<DVect2>()) {
				DVect2 value = vl.at(0).value<DVect2>();
				mingap = value.minComp();
				maxgap = value.maxComp();
			}
			else if (!vl.at(0).isNull())
				throw Exception("gap value %1 not recognized in method bond in contact model %2.", vl.at(0), getName());
			double gap = c->getGap();
			if (gap >= mingap && gap <= maxgap) {
				if (pbProps_)
					pbProps_->pb_state_ = 3;
				else {
					pbProps_ = NEWC(pbProps());
					pbProps_->pb_state_ = 3;
				}
				return true;
			}
			return false;
		}
		case kwPbUnbond: {
			if (!pbProps_ || pbProps_->pb_state_ == 0) return false;
			double mingap = -1.0 * limits<double>::max();
			double maxgap = 0;
			if (vl.at(0).canConvert<double>())
				maxgap = vl.at(0).toDouble();
			else if (vl.at(0).canConvert<DVect2>()) {
				DVect2 value = vl.at(0).value<DVect2>();
				mingap = value.minComp();
				maxgap = value.maxComp();
			}
			else if (!vl.at(0).isNull())
				throw Exception("gap value %1 not recognized in method unbond in contact model %2.", vl.at(0), getName());
			double gap = c->getGap();
			if (gap >= mingap && gap <= maxgap) {
				pbProps_->pb_state_ = 0;
				return true;
			}
			return false;
		}
		case kwArea: {
			if (!userArea_) {
				double rsq(1. / std::max(c->getEnd1Curvature().y(), c->getEnd2Curvature().y()));
#ifdef THREED
				userArea_ = rsq * rsq * dPi;
#else
				userArea_ = rsq * 2.0;
#endif                            
			}
			return true;
		}
		}
		return false;
	}

	double ContactModelmacroelement::getEnergy(uint i) const {
		double ret(0.0);
		if (!energies_)
			return ret;
		switch (i) {
		case kwEStrain:  return energies_->estrain_;
		case kwESlip:    return energies_->eslip_;
		case kwEDashpot: return energies_->edashpot_;
		case kwEPbStrain:return energies_->epbstrain_;
		case kwMeEStrain:  return energies_->me_estrain_;
		case kwMeEPlastic: return energies_->me_eplastic_;
		case kwMeEBStrain: return energies_->me_ebstrain_;
		case kwMeEBPlastic:return energies_->me_ebplastic_;
		}
		assert(0);
		return ret;
	}

	bool ContactModelmacroelement::getEnergyAccumulate(uint i) const {
		switch (i) {
		case kwEStrain:   return false;
		case kwESlip:     return true;
		case kwEDashpot:  return true;
		case kwEPbStrain: return false;
		case kwMeEStrain:  return false;
		case kwMeEPlastic: return true;
		case kwMeEBStrain: return false;
		case kwMeEBPlastic:return true;
		}
		assert(0);
		return false;
	}

	void ContactModelmacroelement::setEnergy(uint i, const double &d) {
		if (!energies_) return;
		switch (i) {
		case kwEStrain:   energies_->estrain_ = d; return;
		case kwESlip:     energies_->eslip_ = d; return;
		case kwEDashpot:  energies_->edashpot_ = d; return;
		case kwEPbStrain: energies_->epbstrain_ = d; return;
		case kwMeEStrain:   energies_->me_estrain_ = d; return;
		case kwMeEPlastic:  energies_->me_eplastic_ = d; return;
		case kwMeEBStrain:  energies_->me_ebstrain_ = d; return;
		case kwMeEBPlastic: energies_->me_ebplastic_ = d; return;
		}
		assert(0);
		return;
	}

	bool ContactModelmacroelement::validate(ContactModelMechanicalState *state, const double &) {
		assert(state);
		const IContactMechanical *c = state->getMechanicalContact();
		assert(c);

		if (state->trackEnergy_)
			activateEnergy();

		if (inheritanceField_ & linKnMask)
			updateKn(c);
		if (inheritanceField_ & linKsMask)
			updateKs(c);
		if (inheritanceField_ & linFricMask)
			updateFric(c);

		updateEffectiveStiffness(state);
		return checkActivity(state->gap_);
	}

	static const QString knstr("kn");
	bool ContactModelmacroelement::updateKn(const IContactMechanical *con) {
		assert(con);
		QVariant v1 = con->getEnd1()->getProperty(knstr);
		QVariant v2 = con->getEnd2()->getProperty(knstr);
		if (!v1.isValid() || !v2.isValid())
			return false;
		double kn1 = v1.toDouble();
		double kn2 = v2.toDouble();
		double val = kn_;
		if (kn1 && kn2)
			kn_ = kn1 * kn2 / (kn1 + kn2);
		else if (kn1)
			kn_ = kn1;
		else if (kn2)
			kn_ = kn2;
		return ((kn_ != val));
	}

	static const QString ksstr("ks");
	bool ContactModelmacroelement::updateKs(const IContactMechanical *con) {
		assert(con);
		QVariant v1 = con->getEnd1()->getProperty(ksstr);
		QVariant v2 = con->getEnd2()->getProperty(ksstr);
		if (!v1.isValid() || !v2.isValid())
			return false;
		double ks1 = v1.toDouble();
		double ks2 = v2.toDouble();
		double val = ks_;
		if (ks1 && ks2)
			ks_ = ks1 * ks2 / (ks1 + ks2);
		else if (ks1)
			ks_ = ks1;
		else if (ks2)
			ks_ = ks2;
		return ((ks_ != val));
	}

	static const QString fricstr("fric");
	bool ContactModelmacroelement::updateFric(const IContactMechanical *con) {
		assert(con);
		QVariant v1 = con->getEnd1()->getProperty(fricstr);
		QVariant v2 = con->getEnd2()->getProperty(fricstr);
		if (!v1.isValid() || !v2.isValid())
			return false;
		double fric1 = std::max(0.0, v1.toDouble());
		double fric2 = std::max(0.0, v2.toDouble());
		double val = fric_;
		fric_ = std::min(fric1, fric2);
		return ((fric_ != val));
	}

	bool ContactModelmacroelement::endPropertyUpdated(const QString &name, const IContactMechanical *c) {
		assert(c);
		QStringList availableProperties = getProperties().simplified().replace(" ", "").split(",", QString::SkipEmptyParts);
		QRegExp rx(name, Qt::CaseInsensitive);
		int idx = availableProperties.indexOf(rx) + 1;
		bool ret = false;

		if (idx <= 0)
			return ret;

		switch (idx) {
		case kwLinKn: { //kn
			if (inheritanceField_ & linKnMask)
				ret = updateKn(c);
			break;
		}
		case kwLinKs: { //ks
			if (inheritanceField_ & linKsMask)
				ret = updateKs(c);
			break;
		}
		case kwLinFric: { //fric
			if (inheritanceField_ & linFricMask)
				updateFric(c);
			break;
		}
		}
		return ret;
	}

	void ContactModelmacroelement::updateEffectiveStiffness(ContactModelMechanicalState *state) {
		DVect2 ret(kn_, ks_);
		// account for viscous damping
		if (dpProps_) {
			DVect2 correct(1.0);
			if (dpProps_->dp_nratio_)
				correct.rx() = sqrt(1.0 + dpProps_->dp_nratio_*dpProps_->dp_nratio_) - dpProps_->dp_nratio_;
			if (dpProps_->dp_sratio_)
				correct.ry() = sqrt(1.0 + dpProps_->dp_sratio_*dpProps_->dp_sratio_) - dpProps_->dp_sratio_;
			ret /= (correct*correct);
		}
		effectiveTranslationalStiffness_ = ret;
		if (isBonded()) {
			double Cmin1 = state->end1Curvature_.x();
			double Cmax1 = state->end1Curvature_.y();
			double Cmax2 = state->end2Curvature_.y();
			double dthick = (Cmin1 == 0.0) ? 1.0 : 0.0;
			double br = pbProps_->pb_rmul_ * 1.0 / std::max(Cmax1, Cmax2);
			if (userArea_)
#ifdef THREED
				br = std::sqrt(userArea_ / dPi);
#else
				br = userArea_ / 2.0;
#endif
			double br2 = br * br;
			double pbArea = dthick <= 0.0 ? dPi * br2 : 2.0*br*dthick;
			double bi = dthick <= 0.0 ? 0.25*pbArea*br2 : 2.0*br*br2*dthick / 3.0;
			pbProps_->pbTransStiff_.rx() = pbProps_->pb_kn_*pbArea;
			pbProps_->pbTransStiff_.ry() = pbProps_->pb_ks_*pbArea;
#if DIM==3 
			pbProps_->pbAngStiff_.rx() = pbProps_->pb_ks_* 2.0 * bi;
			pbProps_->pbAngStiff_.ry() = pbProps_->pb_kn_* bi;
#endif
			pbProps_->pbAngStiff_.rz() = pbProps_->pb_kn_* bi;
		}
	}

	double ContactModelmacroelement::pbStrainEnergy() const {
		double ret(0.0);
		if (pbProps_->pb_kn_)
			ret = 0.5 * pbProps_->pb_F_.x() * pbProps_->pb_F_.x() / pbProps_->pbTransStiff_.x();
		if (pbProps_->pb_ks_) {
			DVect tmp = pbProps_->pb_F_;
			tmp.rx() = 0.0;
			double smag2 = tmp.mag2();
			ret += 0.5 * smag2 / pbProps_->pbTransStiff_.y();
		}

#ifdef THREED
		if (pbProps_->pbAngStiff_.x())
			ret += 0.5 * pbProps_->pb_M_.x() * pbProps_->pb_M_.x() / pbProps_->pbAngStiff_.x();
#endif
		if (pbProps_->pbAngStiff_.z()) {
			DAVect tmp = pbProps_->pb_M_;
#ifdef THREED
			tmp.rx() = 0.0;
			double smag2 = tmp.mag2();
#else
			double smag2 = tmp.z() * tmp.z();
#endif
			ret += 0.5 * smag2 / pbProps_->pbAngStiff_.z();
		}
		return ret;
	}


	double ContactModelmacroelement::calculateElasticModulus(double u_el){

		//return me_ea_;
		// if the wire sliding model is active
		if (me_rata_ > 0 ) {

			double kSlide;
			double kEA = me_ea_;
			if (u_el < me_ratb_)  { // check to avoid extreme values
				kSlide = -me_rata_ / (u_el- me_ratb_);// pow(u_el - me_ratb_, me_rata_);
				kEA = 1.0 / (1.0 / me_ea_ + 1.0 / kSlide);
				me_u_slide(u_el - me_ratb_);
			}

			return kEA;
		}

			return me_ea_;
		
	}

	double ContactModelmacroelement::yieldf(double u_p, double omega_p, double F, double M) {

		double eta = 1.0 - exp(-me_kf_ * abs(u_p) - me_km_ * abs(omega_p));

		me_nf_ = 2.0 * (me_nzero_ + eta * (me_nult_ - me_nzero_));
		double delta_F_ = F - me_nzero_;
		double delta_nf_ = me_nf_ - 2.0 * me_nzero_;
		double F_star = F+me_nf_/2;
		double h = exp(-1.0 / me_betaf_ * pow(-0.5 + F_star / me_nf_,2.0));

		double f;
		double D = me_d_;
		double MD = (M / D); // negative moments fixed
		
		if (F < 0) {
			double MY = D * me_nf_ * sqrt(me_mff_) / 2;
			double MYD = MY / D;
			f = (MYD * (F * F / me_ntens_ / me_ntens_ + MD * MD / MYD / MYD - 1.0));
			double fneg = f;
		}
		else if (F >= 0)
			f = MD * MD + me_mff_ * h * (F_star - me_nf_) * F_star;

		if (isnan(f))
			throw Exception("Nan values encounted while evaluating yield");
		return f;
	}
	// sum the systemStatus properties
	ContactModelmacroelement::systemStatus ContactModelmacroelement::updateStatus(systemStatus myStatus, epIncrements myIncrements) {
	
		myStatus.nf_  += myIncrements.nf_dot_;
		myStatus.u_p_  += myIncrements.u_p_dot_;
		myStatus.omega_p_+= myIncrements.omega_p_dot_;
		myStatus.moment_  += myIncrements.moment_dot_;
		myStatus.force_ += myIncrements.force_dot_;

		return myStatus;
	}

	// multiply the epIncrement properties
	ContactModelmacroelement::epIncrements ContactModelmacroelement::epIncProduct(epIncrements myIncrements, double myMult) {
		myIncrements.nf_dot_*=myMult;
		myIncrements.u_p_dot_*=myMult;
		myIncrements.omega_p_dot_*=myMult;
		myIncrements.moment_dot_*=myMult;
		myIncrements.force_dot_*=myMult;

		return myIncrements;
	}

	// sum the epIncrement properties
	ContactModelmacroelement::epIncrements ContactModelmacroelement::epIncSum(epIncrements myIncrements1, epIncrements myIncrements2) {
		epIncrements newIncrement;
		newIncrement.nf_dot_ = myIncrements1.nf_dot_+myIncrements2.nf_dot_;
		newIncrement.u_p_dot_ = myIncrements1.u_p_dot_ + myIncrements2.u_p_dot_;
		newIncrement.omega_p_dot_ = myIncrements1.omega_p_dot_ + myIncrements2.omega_p_dot_;
		newIncrement.force_dot_ = myIncrements1.force_dot_ + myIncrements2.force_dot_;
		newIncrement.moment_dot_ = myIncrements1.moment_dot_ + myIncrements2.moment_dot_;

		return newIncrement;
	}


	ContactModelmacroelement::systemStatus ContactModelmacroelement::plasticIncrements(systemStatus curr_state, dfIncrements dfinc) {

		// numerical params and vars
		double tiny = 1.0e-13;
		double error_tol = 1.0e-6;
		int ksub = 0;
		int max_ksub = 1000;

		// init
		double T_j = 0.0;
		double DT_j = 1.0;

		systemStatus temp_status = curr_state;
		dfIncrements temp_dfinc = dfinc;

		while (T_j<1.0){


			ksub++;
			
			if (ksub > max_ksub) {
				
				throw Exception("Too many iterations in the plasticity substepping");
				}
			
			// calculate the approximate low order solutions to integrate the elasto-plastic increments using 3rd Order Runge-Kutta
			epIncrements first_approximation = evaluatePlasticity(curr_state, dfinc);
			systemStatus updated_state_2nd_order = updateStatus(curr_state, epIncProduct(first_approximation,DT_j/2.0));

			epIncrements second_approximation = evaluatePlasticity(updated_state_2nd_order, dfinc);
			systemStatus updated_state_3rd_order = updateStatus(curr_state, epIncSum(epIncProduct(first_approximation,-DT_j), epIncProduct(second_approximation,DT_j*2.0)));

			epIncrements third_approximation  = evaluatePlasticity(updated_state_3rd_order, dfinc);

			// form approximate solution of 2nd and 3rd order
			systemStatus state_hat = updateStatus(curr_state, epIncProduct(second_approximation, DT_j));
			systemStatus state_tilde = updateStatus(curr_state, epIncSum(epIncSum(epIncProduct(first_approximation, DT_j * 1.0 / 6.0), epIncProduct(second_approximation, DT_j *2.0 / 3.0)), epIncProduct (third_approximation, DT_j *1.0 / 6.0)));


			// calculate the error residuals
			double sig_tilde[]{state_tilde.force_,state_tilde.moment_}; 
			double nf_tilde = state_tilde.nf_;

			double delta_sig[]{ sig_tilde[0] - state_hat.force_,sig_tilde[1]-state_hat.moment_};
			double delta_nf = nf_tilde - state_hat.nf_;

			double norm_sig = sqrt(std::inner_product(std::begin(sig_tilde), std::end(sig_tilde), std::begin(sig_tilde), 0.0));
			double norm_nf = abs(nf_tilde);


			// assign a minimum value to avoid division by zero errors
			if (norm_sig < tiny)
				norm_sig = tiny;
			if (norm_nf < tiny)
				norm_nf = tiny;

			double vec_residual[3]{ 0.0,0.0,0.0 };
			vec_residual[0] = (1.0 / norm_sig)*delta_sig[0];
			vec_residual[1] = (1.0 / norm_sig)*delta_sig[1];
			vec_residual[2] = (1.0 / norm_nf)*delta_nf;
			double norm_Res = sqrt(std::inner_product(std::begin(vec_residual), std::end(vec_residual), std::begin(vec_residual), 0.0));


			if (norm_Res < tiny)
				norm_Res = tiny;

			// tstep
			double NSS = 0.9 * DT_j*pow(error_tol / norm_Res, 1.0 / 3.0);

			if (norm_Res < error_tol) {
				curr_state = state_tilde;
				T_j = T_j + DT_j;
				DT_j = min(4.0*DT_j, NSS);
				DT_j = min(1.0 - T_j, DT_j);
			}
			else {
				DT_j = max(0.25*DT_j, NSS);
			}

		}
		systemStatus updated_state = curr_state;
		me_nf(updated_state.nf_);
		me_u_p(updated_state.u_p_);
		me_omega_p(updated_state.omega_p_);
		me_force(updated_state.force_);
		me_moment(updated_state.moment_);


		return updated_state;
	}

	ContactModelmacroelement::epIncrements ContactModelmacroelement::evaluatePlasticity(systemStatus myState, dfIncrements dfinc) {
		double tolerance = 1.0e-6;

			// retreive the generalised displacement increments
		double omega_dot = dfinc.omega_dot_;
		double u_dot = dfinc.u_dot_;


		 double f = yieldf(myState.u_p_, myState.omega_p_, myState.force_, myState.moment_);

		 double D = me_d_;


		// wire parameters
		double r = 0.5 * D;
		double L = myState.gap_ + 2.0 * r;
		me_curr_l(L); // update the current contact length

		
		// get the plastic displacements
		double u_p_ = myState.u_p_;
		double omega_p_ = myState.omega_p_;


		// get the system state
		 double F = myState.force_;
		 double M = myState.moment_;
		
		 // retreive the current elastic displacements
		 double u_el_ = (L-me_ref_l_) - u_p_;

		 double kEA = calculateElasticModulus(u_el_);
		 double kEI = me_ei_ * kEA / me_ea_;

		 
		 
		int sign_M = (M >= 0) - (M < 0);
		double MD = (M / D);

		int sign_trial_M = (myState.moment_ + kEI * omega_dot >= 0) - (myState.moment_ + kEI * omega_dot < 0);
		int sign_trial_F = (myState.force_ + kEA * u_dot >= 0) - (myState.force_ + kEA * u_dot < 0);

		double eta = 1.0 - exp(-me_kf_ * abs(myState.u_p_) - me_km_ * abs(myState.omega_p_));


		double nf = 2.0 * (me_nzero_ + eta * (me_nult_ - me_nzero_));
		
		// update the internal variable
		me_nf(nf);

		// compute the partial derivatives
		double dnf_dup = 2.0 * me_kf_ * exp(-me_km_ * abs(myState.omega_p_) - me_kf_ * abs(myState.u_p_)) * (me_nult_ - me_nzero_);
		double dnf_domegap = 2.0 * me_km_ * exp(-me_km_ * abs(myState.omega_p_) - me_kf_ * abs(myState.u_p_)) * (me_nult_ - me_nzero_);

		double F_star = F + nf / 2;

		double df_dM;
		double df_dF;
		double df_dnf;
		double dg_dM;
		double dg_dF;





		int myFsign = (F >= 0) - (F < 0);

		// compressive force domain
		if (F < 0) {

			double MY = nf * D * sqrt(me_mff_) / 2;
			double MYD = MY / D;
			double dmy_dnf = D * sqrt(me_mff_) / 2;
			double tan_theta = MYD / me_ntens_;
			double cos_theta = cos(atan(tan_theta));
			double sin_theta = sin(atan(tan_theta));


			df_dF = 2 * F * MY / (D * me_ntens_ * me_ntens_);
			df_dM = 2 * MD / MY;

			double df_dmy = (F * F / me_ntens_ / me_ntens_ + M * M / MY / MY - 1.0) / D - 2 * M * M / D / MY / MY;

			df_dnf = df_dmy * dmy_dnf;
			dg_dF = df_dF;
			dg_dM = df_dM;

		}
		else // Tensile force domain
		{
			double temp1 = exp(-1.0 / me_betaf_ * pow(F_star / nf - 0.5, 2));
			double temp2 = F - nf / 2;
			
			df_dF = me_mff_* temp1* temp2 + me_mff_ * temp1 * F_star - 1.0 / (nf * me_betaf_) * me_mff_ * temp1 * temp2 * F_star * (F_star / nf - 0.5) * 2;
			df_dM = 2 * M / D / D;
			df_dnf =  me_mff_* temp1* temp2 * 0.5 - me_mff_ * temp1 * F_star * 0.5 + me_mff_ / me_betaf_ * temp1 * temp2 * F_star * (F_star / nf - 0.5) * (F_star / nf / nf - 1.0 / nf * 0.5) * 2;

			dg_dF =  me_mfg_* temp1* temp2 + me_mfg_ * temp1 * F_star - 1.0 / (nf * me_betag_) * me_mfg_ * temp1 * temp2 * F_star * (F_star / nf - me_alphag_) * 2.0;
			dg_dM = 2 * M / D / D;
			
		}

	
		// build the plastic multiplier
		double nomin = df_dM * kEI * omega_dot + df_dF * kEA * u_dot;
		double denom1 = -1.0 * df_dnf * (dnf_domegap * (dg_dM)+dnf_dup * dg_dF); // plastic increment portion
		double denom2 = df_dM * dg_dM * kEI +df_dF * dg_dF * kEA;                // elastic increment portion
		double lambda = nomin / (denom1 + denom2);

		// calculate the plastic increments
		epIncrements myIncrements;

		myIncrements.u_p_dot_ = abs(lambda * dg_dF);
		myIncrements.omega_p_dot_ = abs(lambda * dg_dM);
		myIncrements.force_dot_ = (u_dot - lambda * dg_dF) * kEA;
		myIncrements.moment_dot_ = (omega_dot - lambda * dg_dM) * kEI;

		double nfdot = dnf_dup * abs(myIncrements.u_p_dot_) + dnf_domegap * abs(myIncrements.omega_p_dot_);

		// backcalculating the yield surface size increment (i.e. similar to return mapping, not implemented)
		/*
		double nf2 = 2 * (me_nzero_ - me_nult_ * (exp(-me_kf_ * abs(u_p_ + myIncrements.u_p_dot_) - me_km_ * abs(omega_p_ + myIncrements.omega_p_dot_)) - 1.0));
		double nfdot2 = nf2 - nf;

		eta = 1.0 - exp(-me_kf_ * abs(u_p_ + myIncrements.u_p_dot_) - me_km_ * abs(omega_p_ + myIncrements.omega_p_dot_));
		me_nf_ = 2.0 * (me_nzero_ + eta * (me_nult_ - me_nzero_));
		nf2 = 2.0 * (me_nzero_ + eta * (me_nult_ - me_nzero_));
		double nfdot3 = nf2 - nf;
		*/
		myIncrements.nf_dot_ = nfdot;

		if (isnan(nfdot))
			throw Exception("The yield surface size increment is not a number");

			return myIncrements;


	}
	



		ContactModelmacroelement::systemStatus ContactModelmacroelement::stressStateEvolution(systemStatus initial_state, DVect trans, DAVect ang, double gap, ContactModelMechanicalState *state) {
		double increment_size_tolerance = 1.0e-3;
		double tolerance = 1.0e-3;
		int iteration = 1;
		int max_iter = 1000;

		// wire parameters
		double sa = ((me_d_ / 2) * (me_d_ / 2) * M_PI);
		double im = me_ei_ / (me_ea_ / sa);
		double D = me_d_;
		double im_polar = im * 2.0;
		double beta = me_betaf_;

		double current_length = gap + me_d_;
		me_curr_l(current_length);

		double me_poisson_ = 0.3;

		double curr_force = initial_state.force_;
		double curr_moment = initial_state.moment_;

		double d_u_tot = trans.x();
		double curr_moment_bond = me_moment_ * me_moment_sign_;
		DAVect vec_d_omega_tot = ang*180.0 / dPi;


		double curr_force_D = curr_force;
		
	
		systemStatus curr_state = initial_state;
		systemStatus updated_state;

		
		//////////////////////////////////////////////
		///////////// M CALCULATIONS /////////////////
		//////////////////////////////////////////////
		int sign_curr_moment = me_moment_sign_;

		// the direction of the angle increment (i.e. to have negative bending M)
		int sign_rot_increment = (vec_d_omega_tot.y() >= 0) - (vec_d_omega_tot.y() < 0);
		
		// overall whether it is increasing or decreasing moment
		int sign_yb = sign_curr_moment * sign_rot_increment;

		// at low moment the rotation increment can flip the sign and this is used to keep track of it
		int sign_moment_inversion = 1;

		// then i calculate the modulus of the overall increment and multiply it by the sign calculated just above to obtain the overall increment
		double d_omega_tot = abs(vec_d_omega_tot.y());
		
		// set the initial increment values for the iteration
		double d_u = d_u_tot;
		double d_omega = d_omega_tot;

		 double f;
		 double ftrial;

		 double u_el_ = (current_length - me_ref_l_) - me_u_p_;
		 double kEA = calculateElasticModulus(u_el_);
		 double kEI = me_ei_ * kEA / me_ea_;

		 // use the bisection algorithm to identify the intersection with the yield surface
		 while (abs(d_u) < abs(d_u_tot) || abs(d_omega) < abs(d_omega_tot) || iteration == 1) {

			iteration++;
	
			double trial_force = curr_force + kEA  * d_u;
			double trial_moment = curr_moment + kEI * d_omega * sign_yb;

			if (trial_moment < 0) {
				sign_moment_inversion = -1;
				trial_moment = abs(trial_moment);
			}

			systemStatus trial_state = curr_state;
			trial_state.force_ = trial_force;
			trial_state.moment_ = trial_moment;

			f = yieldf(initial_state.u_p_, initial_state.omega_p_, curr_force, curr_moment);
			ftrial = yieldf(initial_state.u_p_, initial_state.omega_p_, trial_force, trial_moment);


			dfIncrements myIncrements;
			myIncrements.u_dot_ = d_u;
			myIncrements.omega_dot_ = d_omega;

			if (ftrial < 0.0)  //pure elastic increment
				updated_state = trial_state;			
			else if (!me_isplastic() && f < 0.0 && ftrial > 0.0) {
				double a0 = 0.0;
				double a1 = 1.0;
				double a = (a0 + a1) / 2.0;;
				int ksub = 0;
				int maxksub = 1000;


				while (ftrial > 0) {
					a = (a0 + a1) / 2.0;

					trial_force = curr_force + kEA  * d_u * a;
					trial_moment = curr_moment + kEI * d_omega * a;

					ftrial = yieldf(initial_state.u_p_, initial_state.omega_p_, trial_force, trial_moment);

					if (ftrial * f < 0)
						a1 = a;
					else
						a0 = a;

					ksub++;
					if (ksub > maxksub)
						throw Exception("Too many iterations while trying to find the intersection with the yieldf");
				}
			
			if (!me_isplastic()) {
				trial_state.force_ = trial_force;
				trial_state.moment_ = trial_moment;
				curr_force = trial_force;
				curr_moment = trial_moment;
				curr_state = trial_state;
			}
			
				d_u = d_u_tot * (1.0 - a);
				d_omega = d_omega_tot * (1.0 - a);
				myIncrements.omega_dot_ = d_omega;
				myIncrements.u_dot_ = d_u;
				me_isplastic(true);
				iteration = 1;
				updated_state.force_ = trial_force;
				updated_state.moment_ = trial_moment;
			}
			else if (me_isplastic() && ftrial > 0) {
				updated_state = plasticIncrements(curr_state, myIncrements);
				d_u = d_u_tot;
				d_omega = d_omega_tot;
			}
			else {
				char buffer[100];
				sprintf(buffer,"Cannot find the yield surface, f: %.2e, ftrial: $.2e",f,ftrial);
				throw Exception(buffer);
			}
		
			if (iteration > max_iter)
				break;
				
		}


		// i have the updated state in terms of stress increment and bending M increment
		double M_increment = updated_state.moment_-initial_state.moment_;
		double force_increment = updated_state.force_ - initial_state.force_;
		double omegap_increment = updated_state.omega_p_ - initial_state.omega_p_;
		double up_increment = updated_state.u_p_ - initial_state.u_p_;

		if (updated_state.moment_ < 0)
			updated_state.moment_ = 0.0;
			//throw Exception("Negative moment ???");
		if (up_increment<0)
			throw Exception("Negative plastic strain increment?");


		double wp_force = up_increment* (updated_state.force_ + initial_state.force_) / 2;
		if (wp_force < 0)
			wp_force = 0;

		double wp_moment = omegap_increment * (updated_state.moment_ + initial_state.moment_) / 2;
		double delta_wp = wp_force + wp_moment;
		me_wp_ += delta_wp;

		me_wp(me_wp_);

		DVect u_s = trans;
		u_s.rx() = 0;

		double kEA_shear = kEA / (2.0 * (1.0 + me_poisson_));


		DAVect moment_inc;
		moment_inc[0] = ang[0] * 180.0 / dPi * kEA_shear / sa * im_polar;
		moment_inc[1] = M_increment * sign_moment_inversion * sign_curr_moment;
		moment_inc[2] = ang[2] * 180.0/dPi * kEA / sa  * im;

		DVect force_inc;
		force_inc.rx() = force_increment;
		force_inc.ry() = u_s.y()* kEA_shear ;
		force_inc.rz() = u_s.z()* kEA_shear ;

		double damage_multiplier = 1.0;
		double dmgv2 = pow((me_wpl2_ - me_wp_) / (me_wpl2_ - me_wpl1_), me_dalpha_);

		double normalised_wp = (me_wp_ - me_wpl1_) / (me_wpl2_ - me_wpl1_);

		if (me_wp_> me_wpl1_)
			damage_multiplier = 1.0 - pow(normalised_wp, me_dalpha_);
		if (isnan(damage_multiplier)) // if outside the region of interest, set the current force/moment to zero
			damage_multiplier = 0.0;
		if (damage_multiplier < 0.01)
			damage_multiplier = 0.01;
					

		DAVect vec_moment;
		vec_moment[0] = pbProps_->pb_M_[0];
		vec_moment[1] = curr_moment_bond ;
		vec_moment[2] = pbProps_->pb_M_[2];

		DVect vec_force;
		vec_force.rx() = me_force_ ;
		vec_force.ry() = pbProps_->pb_F_[1];
		vec_force.rz() = pbProps_->pb_F_[2];


		pbProps_->pb_F_ -= force_inc;
		pbProps_->pb_M_ -= moment_inc;

		me_force(updated_state.force_);
		me_moment(updated_state.moment_);
		me_force_damaged(updated_state.force_* damage_multiplier);
		me_moment_damaged(updated_state.moment_* damage_multiplier);
		me_moment_sign(sign_curr_moment* sign_moment_inversion);
		me_u_p(updated_state.u_p_);
		me_omega_p(updated_state.omega_p_);
		me_nf(updated_state.nf_);

		return updated_state;

	}

	// this is the function defined by PFC by default, we take its parameters/variables and use them in the macroelement
	bool ContactModelmacroelement::forceDisplacementLaw(ContactModelMechanicalState *state, const double &timestep) {
		assert(state);


		DVect trans = state->relativeTranslationalIncrement_;
		DAVect ang = state->relativeAngularIncrement_;
		double correction = 1.0;
		double overlap = rgap_ - state->gap_;

		if (state->activated()) {
			if (cmEvents_[fActivated] >= 0) {
				auto c = state->getContact();
				std::vector<fish::Parameter> arg = { fish::Parameter(c->getIThing()) };
				IFishCallList* fi = const_cast<IFishCallList*>(state->getProgram()->findInterface<IFishCallList>());
				fi->setCMFishCallArguments(c, arg, cmEvents_[fActivated]);
			}
			if (lin_mode_ == 0 && trans.x()) {
				correction = -1.0 * overlap / trans.x();
				if (correction < 0)
					correction = 1.0;
			}
		}

		systemStatus myStatus;
		myStatus.nf_ = me_nf();
		myStatus.u_p_ = me_u_p();
		myStatus.omega_p_ = me_omega_p();
		myStatus.force_ = me_force();
		myStatus.moment_ = me_moment();
		myStatus.gap_ = state->gap_;

		double previous_force = myStatus.force_;
		double previous_u_p_ = me_u_p();
		double previous_moment = myStatus.moment_;
		double previous_omega_p_ = me_omega_p();

		
		double gap = state->gap_;

		systemStatus updated_state = stressStateEvolution(myStatus, trans, ang,gap,state);
		double damage_multiplier = 1.0;
		double normalised_wp = (me_wp_ - me_wpl1_) / (me_wpl2_ - me_wpl1_);
		if (me_wp_ > me_wpl1_)
			damage_multiplier = 1.0 - pow(normalised_wp, me_dalpha_);
		// mesh failure

		// activate mesh failure
		if (damage_multiplier <= 0.01) {
			pbProps_->pb_state_ = 1;
			pbProps_->pb_F_.fill(0.0);
			pbProps_->pb_M_.fill(0.0);
			if (cmEvents_[fBondBreak] >= 0) {
				auto c = state->getContact();
				double se = -1.0;
				std::vector<fish::Parameter> arg = { fish::Parameter(c->getIThing()),
													 fish::Parameter((qint64)pbProps_->pb_state_),
													 fish::Parameter(pbProps_->pb_ten_),
													 fish::Parameter(se) };
				IFishCallList* fi = const_cast<IFishCallList*>(state->getProgram()->findInterface<IFishCallList>());
				fi->setCMFishCallArguments(c, arg, cmEvents_[fBondBreak]);
			}
		}
		

		if (state->trackEnergy_) {
			assert(energies_);
			energies_->me_estrain_ = -me_force()*me_force() * 0.5 / me_ea_; // maybe need to account for section area and inertia M?
			energies_->me_eplastic_ -= previous_force * me_force() * (me_u_p() - previous_u_p_) * 0.5;
			energies_->me_ebstrain_ = me_moment()*me_moment() * 0.5 / me_ea_;
			energies_->me_ebplastic_ += previous_moment * me_moment() * (me_omega_p() - previous_omega_p_) * 0.5;
		}


		if (state->trackEnergy_) {
			assert(energies_);
		}

		assert(lin_F_ == lin_F_);
		return checkActivity(state->gap_);
	}

	bool ContactModelmacroelement::thermalCoupling(ContactModelMechanicalState *ms, ContactModelThermalState *ts, IContactThermal *ct, const double &) {
		bool ret = false;
		if (!pbProps_) return ret;
		if (pbProps_->pb_state_ < 3) return ret;
		int idx = ct->getModel()->getContactModel()->isProperty("thexp");
		if (idx <= 0) return ret;

		double thexp = (ct->getModel()->getContactModel()->getProperty(idx)).toDouble();
		double length = ts->length_;
		double delTemp = ts->tempInc_;
		double delUn = length * thexp * delTemp;
		if (delUn == 0.0) return ret;

		double dthick = 0.0;
		double Cmin1 = ms->end1Curvature_.x();
		double Cmax1 = ms->end1Curvature_.y();
		double Cmin2 = ms->end2Curvature_.x();
		double Cmax2 = ms->end2Curvature_.y();

		Cmin2;
		if (Cmin1 == 0.0)
			dthick = 1.0;

		double br = pbProps_->pb_rmul_ * 1.0 / std::max(Cmax1, Cmax2);
		if (userArea_)
#ifdef THREED
			br = std::sqrt(userArea_ / dPi);
#else
			br = userArea_ / 2.0;
#endif
		double br2 = br * br;
		double pbArea = dthick <= 0.0 ? dPi * br2 : 2.0*br*dthick;
		
		DVect finc(0.0);
		finc.rx() = 1.0;
		finc *= pbProps_->pb_kn_*pbArea*delUn;
		pbProps_->pb_F_ += finc;

		//ms->force_ += finc;

		// The state force has been updated - update the state with the resulting torques
		//ms->getMechanicalContact()->updateResultingTorquesLocal(ms->force_,&ms->MOn1_,&ms->MOn2_);
		 ret = true;
		 return ret;
		
	}

	void ContactModelmacroelement::setForce(const DVect &v, IContact *c) {
		lin_F(v);
		if (v.x() > 0)
			rgap_ = c->getGap() + v.x() / kn_;
	}

	void ContactModelmacroelement::propagateStateInformation(IContactModelMechanical* old, const CAxes &oldSystem, const CAxes &newSystem) {
		// Only do something if the contact model is of the same type
		if (old->getContactModel()->getName().compare("macroelement", Qt::CaseInsensitive) == 0 && !isBonded()) {
			ContactModelmacroelement *oldCm = (ContactModelmacroelement *)old;
#ifdef THREED
			// Need to rotate just the shear component from oldSystem to newSystem

			// Step 1 - rotate oldSystem so that the normal is the same as the normal of newSystem
			DVect axis = oldSystem.e1() & newSystem.e1();
			double c, ang, s;
			DVect re2;
			if (!checktol(axis.abs().maxComp(), 0.0, 1.0, 1000)) {
				axis = axis.unit();
				c = oldSystem.e1() | newSystem.e1();
				if (c > 0)
					c = std::min(c, 1.0);
				else
					c = std::max(c, -1.0);
				ang = acos(c);
				s = sin(ang);
				double t = 1. - c;
				DMatrix<3, 3> rm;
				rm.get(0, 0) = t * axis.x()*axis.x() + c;
				rm.get(0, 1) = t * axis.x()*axis.y() - axis.z()*s;
				rm.get(0, 2) = t * axis.x()*axis.z() + axis.y()*s;
				rm.get(1, 0) = t * axis.x()*axis.y() + axis.z()*s;
				rm.get(1, 1) = t * axis.y()*axis.y() + c;
				rm.get(1, 2) = t * axis.y()*axis.z() - axis.x()*s;
				rm.get(2, 0) = t * axis.x()*axis.z() - axis.y()*s;
				rm.get(2, 1) = t * axis.y()*axis.z() + axis.x()*s;
				rm.get(2, 2) = t * axis.z()*axis.z() + c;
				re2 = rm * oldSystem.e2();
			}
			else
				re2 = oldSystem.e2();
			// Step 2 - get the angle between the oldSystem rotated shear and newSystem shear
			axis = re2 & newSystem.e2();
			DVect2 tpf;
			DMatrix<2, 2> m;
			if (!checktol(axis.abs().maxComp(), 0.0, 1.0, 1000)) {
				axis = axis.unit();
				c = re2 | newSystem.e2();
				if (c > 0)
					c = std::min(c, 1.0);
				else
					c = std::max(c, -1.0);
				ang = acos(c);
				if (!checktol(axis.x(), newSystem.e1().x(), 1.0, 100))
					ang *= -1;
				s = sin(ang);
				m.get(0, 0) = c;
				m.get(1, 0) = s;
				m.get(0, 1) = -m.get(1, 0);
				m.get(1, 1) = m.get(0, 0);
				tpf = m * DVect2(oldCm->lin_F_.y(), oldCm->lin_F_.z());
			}
			else {
				m.get(0, 0) = 1.;
				m.get(0, 1) = 0.;
				m.get(1, 0) = 0.;
				m.get(1, 1) = 1.;
				tpf = DVect2(oldCm->lin_F_.y(), oldCm->lin_F_.z());
			}
			DVect pforce = DVect(0, tpf.x(), tpf.y());
#else
			oldSystem;
			newSystem;
			DVect pforce = DVect(0, oldCm->lin_F_.y());
#endif
			for (int i = 1; i < dim; ++i)
				lin_F_.rdof(i) += pforce.dof(i);
			oldCm->lin_F_ = DVect(0.0);
			if (dpProps_ && oldCm->dpProps_) {
#ifdef THREED
				tpf = m * DVect2(oldCm->dpProps_->dp_F_.y(), oldCm->dpProps_->dp_F_.z());
				pforce = DVect(oldCm->dpProps_->dp_F_.x(), tpf.x(), tpf.y());
#else
				pforce = oldCm->dpProps_->dp_F_;
#endif
				dpProps_->dp_F_ += pforce;
				oldCm->dpProps_->dp_F_ = DVect(0.0);
			}
			if (oldCm->getEnergyActivated()) {
				activateEnergy();
				energies_->estrain_ = oldCm->energies_->estrain_;
				energies_->eslip_ = oldCm->energies_->eslip_;
				energies_->edashpot_ = oldCm->energies_->edashpot_;
				energies_->epbstrain_ = oldCm->energies_->epbstrain_;
				oldCm->energies_->estrain_ = 0.0;
				oldCm->energies_->edashpot_ = 0.0;
				oldCm->energies_->eslip_ = 0.0;
				oldCm->energies_->epbstrain_ = 0.0;
			}
			rgap_ = oldCm->rgap_;
		}
		assert(lin_F_ == lin_F_);
	}

	void ContactModelmacroelement::setNonForcePropsFrom(IContactModel *old) {
		// Only do something if the contact model is of the same type
		if (old->getName().compare("macroelement", Qt::CaseInsensitive) == 0 && !isBonded()) {
			ContactModelmacroelement *oldCm = (ContactModelmacroelement *)old;
			kn_ = oldCm->kn_;
			ks_ = oldCm->ks_;
			fric_ = oldCm->fric_;
			lin_mode_ = oldCm->lin_mode_;
			rgap_ = oldCm->rgap_;
			userArea_ = oldCm->userArea_;

			if (oldCm->dpProps_) {
				if (!dpProps_)
					dpProps_ = NEWC(dpProps());
				dpProps_->dp_nratio_ = oldCm->dpProps_->dp_nratio_;
				dpProps_->dp_sratio_ = oldCm->dpProps_->dp_sratio_;
				dpProps_->dp_mode_ = oldCm->dpProps_->dp_mode_;
			}

			if (oldCm->pbProps_) {
				if (!pbProps_)
					pbProps_ = NEWC(pbProps());
				pbProps_->pb_rmul_ = oldCm->pbProps_->pb_rmul_;
				pbProps_->pb_kn_ = oldCm->pbProps_->pb_kn_;
				pbProps_->pb_ks_ = oldCm->pbProps_->pb_ks_;
				pbProps_->pb_mcf_ = oldCm->pbProps_->pb_mcf_;
				pbProps_->pb_fa_ = oldCm->pbProps_->pb_fa_;
				pbProps_->pb_state_ = oldCm->pbProps_->pb_state_;
				pbProps_->pb_coh_ = oldCm->pbProps_->pb_coh_;
				pbProps_->pb_ten_ = oldCm->pbProps_->pb_ten_;
				pbProps_->pbTransStiff_ = oldCm->pbProps_->pbTransStiff_;
				pbProps_->pbAngStiff_ = oldCm->pbProps_->pbAngStiff_;
			}
		}
	}

	DVect ContactModelmacroelement::getForce(const IContactMechanical *) const {
		DVect ret(lin_F_);
		if (dpProps_)
			ret += dpProps_->dp_F_;
		if (pbProps_)
			ret += pbProps_->pb_F_;
		return ret;
	}

	DAVect ContactModelmacroelement::getMOn1(const IContactMechanical *c) const {
		DVect force = getForce(c);
		DAVect ret(0.0);
		if (pbProps_)
			ret = pbProps_->pb_M_;
		c->updateResultingTorqueOn1Local(force, &ret);
		return ret;
	}

	DAVect ContactModelmacroelement::getMOn2(const IContactMechanical *c) const {
		DVect force = getForce(c);
		DAVect ret(0.0);
		if (pbProps_)
			ret = pbProps_->pb_M_;
		c->updateResultingTorqueOn2Local(force, &ret);
		return ret;
	}

	DAVect ContactModelmacroelement::getMomentOn1(const IContactMechanical* c) const {
		DVect force = getForce(c);
	//	force.rx() = 0.0;
	//	force.ry() = 0.0;
	//	force.rz() = 1.0;
		DAVect ret(0.0);
		if (pbProps_)
			ret = pbProps_->pb_M_;
		c->updateResultingTorqueOn1Local(force, &ret);
		return ret;
	}

	DAVect ContactModelmacroelement::getMomentOn2(const IContactMechanical* c) const {
		DVect force = getForce(c);
	//	force.rx() = 0.0;
	//	force.ry() = 0.0;
	//	force.rz() = 1.0;
		DAVect ret(0.0);
		if (pbProps_)
			ret = pbProps_->pb_M_;
		c->updateResultingTorqueOn2Local(force, &ret);
		return ret;
	}

	DVect3 ContactModelmacroelement::pbData(const IContactMechanical *c) const {
		double Cmin1 = c->getEnd1Curvature().x();
		double Cmax1 = c->getEnd1Curvature().y();
		double Cmax2 = c->getEnd2Curvature().y();
		double dthick = (Cmin1 == 0.0) ? c->getEnd1Extent().x() : 0.0;
		double br = pbProps_->pb_rmul_ * 1.0 / std::max(Cmax1, Cmax2);
		if (userArea_)
#ifdef THREED
			br = std::sqrt(userArea_ / dPi);
#else
			br = userArea_ / 2.0;
#endif
		double br2 = br * br;
		double pbArea = dthick <= 0.0 ? dPi * br2 : 2.0*br*dthick;
		double bi = dthick <= 0.0 ? 0.25*pbArea*br2 : 2.0*br*br2*dthick / 3.0;
		return DVect3(pbArea, bi, br);
	}

	DVect2 ContactModelmacroelement::pbSMax(const IContactMechanical *c) const {
		DVect3 data = pbData(c);
		double pbArea = data.x();
		double bi = data.y();
		double br = data.z();
		/* maximum stresses */
		double dbend = sqrt(pbProps_->pb_M_.y()*pbProps_->pb_M_.y() + pbProps_->pb_M_.z()*pbProps_->pb_M_.z());
		double dtwist = pbProps_->pb_M_.x();
		DVect bfs(pbProps_->pb_F_);
		bfs.rx() = 0.0;
		double dbfs = bfs.mag();
		double nsmax = -(pbProps_->pb_F_.x() / pbArea) + pbProps_->pb_mcf_ * dbend * br / bi;
		double ssmax = dbfs / pbArea + pbProps_->pb_mcf_ * std::abs(dtwist) * 0.5* br / bi;
		return DVect2(nsmax, ssmax);
	}

	double ContactModelmacroelement::pbShearStrength(const double &pbArea) const {
		if (!pbProps_) return 0.0;
		double sig = -1.0*pbProps_->pb_F_.x() / pbArea;
		double nstr = pbProps_->pb_state_ > 2 ? pbProps_->pb_ten_ : 0.0;
		return sig <= nstr ? pbProps_->pb_coh_ - std::tan(dDegrad*pbProps_->pb_fa_)*sig
			: pbProps_->pb_coh_ - std::tan(dDegrad*pbProps_->pb_fa_)*nstr;
	}

	void ContactModelmacroelement::setDampCoefficients(const double &mass, double *vcn, double *vcs) {
		*vcn = dpProps_->dp_nratio_ * 2.0 * sqrt(mass*(kn_));
		*vcs = dpProps_->dp_sratio_ * 2.0 * sqrt(mass*(ks_));
	}

} // namespace itascaxd
// EoF