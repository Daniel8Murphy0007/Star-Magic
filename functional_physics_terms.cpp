// ============================================
// FUNCTIONAL PhysicsTerm Wrapper Classes
// Generated: 2025-11-13 04:02:40
// Each wrapper calls its embedded module's compute method
// Total: 126 wrapper classes
// ============================================

// Wrapper for MagnetarSGR0501_4516 (Source14)
class MagnetarSGR0501_4516Term : public PhysicsTerm {
private:
    MagnetarSGR0501_4516 instance;
public:
    MagnetarSGR0501_4516Term() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.compute_Ug(t);
    }
    
    std::string getName() const override {
        return "MagnetarSGR0501_4516";
    }
    
    std::string getDescription() const override {
        return "Source14: MagnetarSGR0501_4516.compute_Ug()";
    }
};

// Wrapper for SMBHSgrAStar (Source15)
class SMBHSgrAStarTerm : public PhysicsTerm {
private:
    SMBHSgrAStar instance;
public:
    SMBHSgrAStarTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.compute_Ug(t);
    }
    
    std::string getName() const override {
        return "SMBHSgrAStar";
    }
    
    std::string getDescription() const override {
        return "Source15: SMBHSgrAStar.compute_Ug()";
    }
};

// Wrapper for AndromedaUQFFModule (Source28)
class AndromedaUQFFModuleTerm : public PhysicsTerm {
private:
    AndromedaUQFFModule instance;
public:
    AndromedaUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "AndromedaUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source28: AndromedaUQFFModule.computeG()";
    }
};

// Wrapper for SombreroUQFFModule (Source29)
class SombreroUQFFModuleTerm : public PhysicsTerm {
private:
    SombreroUQFFModule instance;
public:
    SombreroUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "SombreroUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source29: SombreroUQFFModule.computeG()";
    }
};

// Wrapper for SaturnUQFFModule (Source30)
class SaturnUQFFModuleTerm : public PhysicsTerm {
private:
    SaturnUQFFModule instance;
public:
    SaturnUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "SaturnUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source30: SaturnUQFFModule.computeG()";
    }
};

// Wrapper for M16UQFFModule (Source31)
class M16UQFFModuleTerm : public PhysicsTerm {
private:
    M16UQFFModule instance;
public:
    M16UQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "M16UQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source31: M16UQFFModule.computeG()";
    }
};

// Wrapper for CrabUQFFModule (Source32)
class CrabUQFFModuleTerm : public PhysicsTerm {
private:
    CrabUQFFModule instance;
public:
    CrabUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "CrabUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source32: CrabUQFFModule.computeG()";
    }
};

// Wrapper for SGR1745UQFFModule (Source33)
class SGR1745UQFFModuleTerm : public PhysicsTerm {
private:
    SGR1745UQFFModule instance;
public:
    SGR1745UQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "SGR1745UQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source33: SGR1745UQFFModule.computeG()";
    }
};

// Wrapper for SGR1745UQFFModule (Source34)
class SGR1745UQFFModuleTerm : public PhysicsTerm {
private:
    SGR1745UQFFModule instance;
public:
    SGR1745UQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "SGR1745UQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source34: SGR1745UQFFModule.computeG()";
    }
};

// Wrapper for SgrA_UQFFModule (Source35)
class SgrA_UQFFModuleTerm : public PhysicsTerm {
private:
    SgrA_UQFFModule instance;
public:
    SgrA_UQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "SgrA_UQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source35: SgrA_UQFFModule.computeG()";
    }
};

// Wrapper for TapestryUQFFModule (Source36)
class TapestryUQFFModuleTerm : public PhysicsTerm {
private:
    TapestryUQFFModule instance;
public:
    TapestryUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "TapestryUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source36: TapestryUQFFModule.computeG()";
    }
};

// Wrapper for ResonanceSuperconductiveUQFFModule (Source37)
class ResonanceSuperconductiveUQFFModuleTerm : public PhysicsTerm {
private:
    ResonanceSuperconductiveUQFFModule instance;
public:
    ResonanceSuperconductiveUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeResonanceTerm(t);
    }
    
    std::string getName() const override {
        return "ResonanceSuperconductiveUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source37: ResonanceSuperconductiveUQFFModule.computeResonanceTerm()";
    }
};

// Wrapper for CompressedResonanceUQFFModule (Source38)
class CompressedResonanceUQFFModuleTerm : public PhysicsTerm {
private:
    CompressedResonanceUQFFModule instance;
public:
    CompressedResonanceUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeCompressedResTerm(t);
    }
    
    std::string getName() const override {
        return "CompressedResonanceUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source38: CompressedResonanceUQFFModule.computeCompressedResTerm()";
    }
};

// Wrapper for CrabResonanceUQFFModule (Source39)
class CrabResonanceUQFFModuleTerm : public PhysicsTerm {
private:
    CrabResonanceUQFFModule instance;
public:
    CrabResonanceUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "CrabResonanceUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source39: CrabResonanceUQFFModule.computeG()";
    }
};

// Wrapper for CompressedResonanceUQFF24Module (Source40)
class CompressedResonanceUQFF24ModuleTerm : public PhysicsTerm {
private:
    CompressedResonanceUQFF24Module instance;
public:
    CompressedResonanceUQFF24ModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeCompressedResTerm(t);
    }
    
    std::string getName() const override {
        return "CompressedResonanceUQFF24Module";
    }
    
    std::string getDescription() const override {
        return "Source40: CompressedResonanceUQFF24Module.computeCompressedResTerm()";
    }
};

// Wrapper for UniverseDiameterUQFFModule (Source41)
class UniverseDiameterUQFFModuleTerm : public PhysicsTerm {
private:
    UniverseDiameterUQFFModule instance;
public:
    UniverseDiameterUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "UniverseDiameterUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source41: UniverseDiameterUQFFModule.computeG()";
    }
};

// Wrapper for HydrogenAtomUQFFModule (Source42)
class HydrogenAtomUQFFModuleTerm : public PhysicsTerm {
private:
    HydrogenAtomUQFFModule instance;
public:
    HydrogenAtomUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "HydrogenAtomUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source42: HydrogenAtomUQFFModule.computeG()";
    }
};

// Wrapper for HydrogenPToEResonanceUQFFModule (Source43)
class HydrogenPToEResonanceUQFFModuleTerm : public PhysicsTerm {
private:
    HydrogenPToEResonanceUQFFModule instance;
public:
    HydrogenPToEResonanceUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeResonanceTerm(t);
    }
    
    std::string getName() const override {
        return "HydrogenPToEResonanceUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source43: HydrogenPToEResonanceUQFFModule.computeResonanceTerm()";
    }
};

// Wrapper for LagoonUQFFModule (Source44)
class LagoonUQFFModuleTerm : public PhysicsTerm {
private:
    LagoonUQFFModule instance;
public:
    LagoonUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "LagoonUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source44: LagoonUQFFModule.computeG()";
    }
};

// Wrapper for SpiralSupernovaeUQFFModule (Source45)
class SpiralSupernovaeUQFFModuleTerm : public PhysicsTerm {
private:
    SpiralSupernovaeUQFFModule instance;
public:
    SpiralSupernovaeUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "SpiralSupernovaeUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source45: SpiralSupernovaeUQFFModule.computeG()";
    }
};

// Wrapper for NGC6302UQFFModule (Source46)
class NGC6302UQFFModuleTerm : public PhysicsTerm {
private:
    NGC6302UQFFModule instance;
public:
    NGC6302UQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "NGC6302UQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source46: NGC6302UQFFModule.computeG()";
    }
};

// Wrapper for NGC6302ResonanceUQFFModule (Source47)
class NGC6302ResonanceUQFFModuleTerm : public PhysicsTerm {
private:
    NGC6302ResonanceUQFFModule instance;
public:
    NGC6302ResonanceUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "NGC6302ResonanceUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source47: NGC6302ResonanceUQFFModule.computeG()";
    }
};

// Wrapper for OrionUQFFModule (Source48)
class OrionUQFFModuleTerm : public PhysicsTerm {
private:
    OrionUQFFModule instance;
public:
    OrionUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "OrionUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source48: OrionUQFFModule.computeG()";
    }
};

// Wrapper for CompressedResonanceUQFF34Module (Source49)
class CompressedResonanceUQFF34ModuleTerm : public PhysicsTerm {
private:
    CompressedResonanceUQFF34Module instance;
public:
    CompressedResonanceUQFF34ModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeCompressed(t);
    }
    
    std::string getName() const override {
        return "CompressedResonanceUQFF34Module";
    }
    
    std::string getDescription() const override {
        return "Source49: CompressedResonanceUQFF34Module.computeCompressed()";
    }
};

// Wrapper for MultiUQFFModule (Source52)
class MultiUQFFModuleTerm : public PhysicsTerm {
private:
    MultiUQFFModule instance;
public:
    MultiUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "MultiUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source52: MultiUQFFModule.computeG()";
    }
};

// Wrapper for YoungStarsOutflowsUQFFModule (Source54)
class YoungStarsOutflowsUQFFModuleTerm : public PhysicsTerm {
private:
    YoungStarsOutflowsUQFFModule instance;
public:
    YoungStarsOutflowsUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "YoungStarsOutflowsUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source54: YoungStarsOutflowsUQFFModule.computeG()";
    }
};

// Wrapper for BigBangGravityUQFFModule (Source56)
class BigBangGravityUQFFModuleTerm : public PhysicsTerm {
private:
    BigBangGravityUQFFModule instance;
public:
    BigBangGravityUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "BigBangGravityUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source56: BigBangGravityUQFFModule.computeG()";
    }
};

// Wrapper for MultiCompressedUQFFModule (Source57)
class MultiCompressedUQFFModuleTerm : public PhysicsTerm {
private:
    MultiCompressedUQFFModule instance;
public:
    MultiCompressedUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "MultiCompressedUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source57: MultiCompressedUQFFModule.computeG()";
    }
};

// Wrapper for MultiUQFFCompressionModule (Source60)
class MultiUQFFCompressionModuleTerm : public PhysicsTerm {
private:
    MultiUQFFCompressionModule instance;
public:
    MultiUQFFCompressionModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "MultiUQFFCompressionModule";
    }
    
    std::string getDescription() const override {
        return "Source60: MultiUQFFCompressionModule.computeG()";
    }
};

// Wrapper for UFEOrbModule (Source64)
class UFEOrbModuleTerm : public PhysicsTerm {
private:
    UFEOrbModule instance;
public:
    UFEOrbModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeUP(t);
    }
    
    std::string getName() const override {
        return "UFEOrbModule";
    }
    
    std::string getDescription() const override {
        return "Source64: UFEOrbModule.computeUP()";
    }
};

// Wrapper for NebularUQFFModule (Source65)
class NebularUQFFModuleTerm : public PhysicsTerm {
private:
    NebularUQFFModule instance;
public:
    NebularUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeElectricField(t);
    }
    
    std::string getName() const override {
        return "NebularUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source65: NebularUQFFModule.computeElectricField()";
    }
};

// Wrapper for RedDwarfUQFFModule (Source66)
class RedDwarfUQFFModuleTerm : public PhysicsTerm {
private:
    RedDwarfUQFFModule instance;
public:
    RedDwarfUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeWmag(t);
    }
    
    std::string getName() const override {
        return "RedDwarfUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source66: RedDwarfUQFFModule.computeWmag()";
    }
};

// Wrapper for InertiaUQFFModule (Source67)
class InertiaUQFFModuleTerm : public PhysicsTerm {
private:
    InertiaUQFFModule instance;
public:
    InertiaUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeTwistPhase(t);
    }
    
    std::string getName() const override {
        return "InertiaUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source67: InertiaUQFFModule.computeTwistPhase()";
    }
};

// Wrapper for HydrogenUQFFModule (Source68)
class HydrogenUQFFModuleTerm : public PhysicsTerm {
private:
    HydrogenUQFFModule instance;
public:
    HydrogenUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeEspace(t);
    }
    
    std::string getName() const override {
        return "HydrogenUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source68: HydrogenUQFFModule.computeEspace()";
    }
};

// Wrapper for UQFFCompressionModule (Source69)
class UQFFCompressionModuleTerm : public PhysicsTerm {
private:
    UQFFCompressionModule instance;
public:
    UQFFCompressionModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "UQFFCompressionModule";
    }
    
    std::string getDescription() const override {
        return "Source69: UQFFCompressionModule.computeG()";
    }
};

// Wrapper for M51UQFFModule (Source70)
class M51UQFFModuleTerm : public PhysicsTerm {
private:
    M51UQFFModule instance;
public:
    M51UQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "M51UQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source70: M51UQFFModule.computeG()";
    }
};

// Wrapper for NGC1316UQFFModule (Source71)
class NGC1316UQFFModuleTerm : public PhysicsTerm {
private:
    NGC1316UQFFModule instance;
public:
    NGC1316UQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "NGC1316UQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source71: NGC1316UQFFModule.computeG()";
    }
};

// Wrapper for V838MonUQFFModule (Source72)
class V838MonUQFFModuleTerm : public PhysicsTerm {
private:
    V838MonUQFFModule instance;
public:
    V838MonUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeIecho(t);
    }
    
    std::string getName() const override {
        return "V838MonUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source72: V838MonUQFFModule.computeIecho()";
    }
};

// Wrapper for NGC1300UQFFModule (Source73)
class NGC1300UQFFModuleTerm : public PhysicsTerm {
private:
    NGC1300UQFFModule instance;
public:
    NGC1300UQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "NGC1300UQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source73: NGC1300UQFFModule.computeG()";
    }
};

// Wrapper for UQFFCompressedResonanceModule (Source74)
class UQFFCompressedResonanceModuleTerm : public PhysicsTerm {
private:
    UQFFCompressedResonanceModule instance;
public:
    UQFFCompressedResonanceModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "UQFFCompressedResonanceModule";
    }
    
    std::string getDescription() const override {
        return "Source74: UQFFCompressedResonanceModule.computeG()";
    }
};

// Wrapper for NGC2264UQFFModule (Source76)
class NGC2264UQFFModuleTerm : public PhysicsTerm {
private:
    NGC2264UQFFModule instance;
public:
    NGC2264UQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "NGC2264UQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source76: NGC2264UQFFModule.computeG()";
    }
};

// Wrapper for UGC10214UQFFModule (Source77)
class UGC10214UQFFModuleTerm : public PhysicsTerm {
private:
    UGC10214UQFFModule instance;
public:
    UGC10214UQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "UGC10214UQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source77: UGC10214UQFFModule.computeG()";
    }
};

// Wrapper for NGC4676UQFFModule (Source78)
class NGC4676UQFFModuleTerm : public PhysicsTerm {
private:
    NGC4676UQFFModule instance;
public:
    NGC4676UQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "NGC4676UQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source78: NGC4676UQFFModule.computeG()";
    }
};

// Wrapper for RedSpiderUQFFModule (Source79)
class RedSpiderUQFFModuleTerm : public PhysicsTerm {
private:
    RedSpiderUQFFModule instance;
public:
    RedSpiderUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "RedSpiderUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source79: RedSpiderUQFFModule.computeG()";
    }
};

// Wrapper for SMBHBinaryUQFFModule (Source80)
class SMBHBinaryUQFFModuleTerm : public PhysicsTerm {
private:
    SMBHBinaryUQFFModule instance;
public:
    SMBHBinaryUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "SMBHBinaryUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source80: SMBHBinaryUQFFModule.computeG()";
    }
};

// Wrapper for NGC346UQFFModule (Source81)
class NGC346UQFFModuleTerm : public PhysicsTerm {
private:
    NGC346UQFFModule instance;
public:
    NGC346UQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "NGC346UQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source81: NGC346UQFFModule.computeG()";
    }
};

// Wrapper for SMBHUQFFModule (Source82)
class SMBHUQFFModuleTerm : public PhysicsTerm {
private:
    SMBHUQFFModule instance;
public:
    SMBHUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "SMBHUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source82: SMBHUQFFModule.computeG()";
    }
};

// Wrapper for LENRUQFFModule (Source83)
class LENRUQFFModuleTerm : public PhysicsTerm {
private:
    LENRUQFFModule instance;
public:
    LENRUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeNeutronRate(t);
    }
    
    std::string getName() const override {
        return "LENRUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source83: LENRUQFFModule.computeNeutronRate()";
    }
};

// Wrapper for LENRCalibUQFFModule (Source84)
class LENRCalibUQFFModuleTerm : public PhysicsTerm {
private:
    LENRCalibUQFFModule instance;
public:
    LENRCalibUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeEta(t);
    }
    
    std::string getName() const override {
        return "LENRCalibUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source84: LENRCalibUQFFModule.computeEta()";
    }
};

// Wrapper for NGC346UQFFModule (Source85)
class NGC346UQFFModuleTerm : public PhysicsTerm {
private:
    NGC346UQFFModule instance;
public:
    NGC346UQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "NGC346UQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source85: NGC346UQFFModule.computeG()";
    }
};

// Wrapper for MUGEModule (Source86)
class MUGEModuleTerm : public PhysicsTerm {
private:
    MUGEModule instance;
public:
    MUGEModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG_compressed(t);
    }
    
    std::string getName() const override {
        return "MUGEModule";
    }
    
    std::string getDescription() const override {
        return "Source86: MUGEModule.computeG_compressed()";
    }
};

// Wrapper for MUGEResonanceModule (Source87)
class MUGEResonanceModuleTerm : public PhysicsTerm {
private:
    MUGEResonanceModule instance;
public:
    MUGEResonanceModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG_resonance(t);
    }
    
    std::string getName() const override {
        return "MUGEResonanceModule";
    }
    
    std::string getDescription() const override {
        return "Source87: MUGEResonanceModule.computeG_resonance()";
    }
};

// Wrapper for AndromedaUQFFModule (Source88)
class AndromedaUQFFModuleTerm : public PhysicsTerm {
private:
    AndromedaUQFFModule instance;
public:
    AndromedaUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeG(t);
    }
    
    std::string getName() const override {
        return "AndromedaUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source88: AndromedaUQFFModule.computeG()";
    }
};

// Wrapper for AetherCouplingModule (Source89)
class AetherCouplingModuleTerm : public PhysicsTerm {
private:
    AetherCouplingModule instance;
public:
    AetherCouplingModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeT_s(t);
    }
    
    std::string getName() const override {
        return "AetherCouplingModule";
    }
    
    std::string getDescription() const override {
        return "Source89: AetherCouplingModule.computeT_s()";
    }
};

// Wrapper for BackgroundAetherModule (Source90)
class BackgroundAetherModuleTerm : public PhysicsTerm {
private:
    BackgroundAetherModule instance;
public:
    BackgroundAetherModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeT_s(t);
    }
    
    std::string getName() const override {
        return "BackgroundAetherModule";
    }
    
    std::string getDescription() const override {
        return "Source90: BackgroundAetherModule.computeT_s()";
    }
};

// Wrapper for DPMModule (Source91)
class DPMModuleTerm : public PhysicsTerm {
private:
    DPMModule instance;
public:
    DPMModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeSCmEnergy(t);
    }
    
    std::string getName() const override {
        return "DPMModule";
    }
    
    std::string getDescription() const override {
        return "Source91: DPMModule.computeSCmEnergy()";
    }
};

// Wrapper for BuoyancyCouplingModule (Source92)
class BuoyancyCouplingModuleTerm : public PhysicsTerm {
private:
    BuoyancyCouplingModule instance;
public:
    BuoyancyCouplingModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeBeta(t);
    }
    
    std::string getName() const override {
        return "BuoyancyCouplingModule";
    }
    
    std::string getDescription() const override {
        return "Source92: BuoyancyCouplingModule.computeBeta()";
    }
};

// Wrapper for SolarWindBuoyancyModule (Source93)
class SolarWindBuoyancyModuleTerm : public PhysicsTerm {
private:
    SolarWindBuoyancyModule instance;
public:
    SolarWindBuoyancyModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeEpsilon_sw(t);
    }
    
    std::string getName() const override {
        return "SolarWindBuoyancyModule";
    }
    
    std::string getDescription() const override {
        return "Source93: SolarWindBuoyancyModule.computeEpsilon_sw()";
    }
};

// Wrapper for UgCouplingModule (Source94)
class UgCouplingModuleTerm : public PhysicsTerm {
private:
    UgCouplingModule instance;
public:
    UgCouplingModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeK_i(t);
    }
    
    std::string getName() const override {
        return "UgCouplingModule";
    }
    
    std::string getDescription() const override {
        return "Source94: UgCouplingModule.computeK_i()";
    }
};

// Wrapper for MagneticStringModule (Source95)
class MagneticStringModuleTerm : public PhysicsTerm {
private:
    MagneticStringModule instance;
public:
    MagneticStringModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeRj(t);
    }
    
    std::string getName() const override {
        return "MagneticStringModule";
    }
    
    std::string getDescription() const override {
        return "Source95: MagneticStringModule.computeRj()";
    }
};

// Wrapper for GalacticDistanceModule (Source96)
class GalacticDistanceModuleTerm : public PhysicsTerm {
private:
    GalacticDistanceModule instance;
public:
    GalacticDistanceModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeDg(t);
    }
    
    std::string getName() const override {
        return "GalacticDistanceModule";
    }
    
    std::string getDescription() const override {
        return "Source96: GalacticDistanceModule.computeDg()";
    }
};

// Wrapper for FeedbackFactorModule (Source97)
class FeedbackFactorModuleTerm : public PhysicsTerm {
private:
    FeedbackFactorModule instance;
public:
    FeedbackFactorModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF_feedback(t);
    }
    
    std::string getName() const override {
        return "FeedbackFactorModule";
    }
    
    std::string getDescription() const override {
        return "Source97: FeedbackFactorModule.computeF_feedback()";
    }
};

// Wrapper for UnifiedFieldModule (Source98)
class UnifiedFieldModuleTerm : public PhysicsTerm {
private:
    UnifiedFieldModule instance;
public:
    UnifiedFieldModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeFU(t);
    }
    
    std::string getName() const override {
        return "UnifiedFieldModule";
    }
    
    std::string getDescription() const override {
        return "Source98: UnifiedFieldModule.computeFU()";
    }
};

// Wrapper for HeavisideFractionModule (Source100)
class HeavisideFractionModuleTerm : public PhysicsTerm {
private:
    HeavisideFractionModule instance;
public:
    HeavisideFractionModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF_Heaviside(t);
    }
    
    std::string getName() const override {
        return "HeavisideFractionModule";
    }
    
    std::string getDescription() const override {
        return "Source100: HeavisideFractionModule.computeF_Heaviside()";
    }
};

// Wrapper for HeliosphereThicknessModule (Source101)
class HeliosphereThicknessModuleTerm : public PhysicsTerm {
private:
    HeliosphereThicknessModule instance;
public:
    HeliosphereThicknessModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeH_SCm(t);
    }
    
    std::string getName() const override {
        return "HeliosphereThicknessModule";
    }
    
    std::string getDescription() const override {
        return "Source101: HeliosphereThicknessModule.computeH_SCm()";
    }
};

// Wrapper for UgIndexModule (Source102)
class UgIndexModuleTerm : public PhysicsTerm {
private:
    UgIndexModule instance;
public:
    UgIndexModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeU_gi(t);
    }
    
    std::string getName() const override {
        return "UgIndexModule";
    }
    
    std::string getDescription() const override {
        return "Source102: UgIndexModule.computeU_gi()";
    }
};

// Wrapper for InertiaCouplingModule (Source103)
class InertiaCouplingModuleTerm : public PhysicsTerm {
private:
    InertiaCouplingModule instance;
public:
    InertiaCouplingModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeLambda_i(t);
    }
    
    std::string getName() const override {
        return "InertiaCouplingModule";
    }
    
    std::string getDescription() const override {
        return "Source103: InertiaCouplingModule.computeLambda_i()";
    }
};

// Wrapper for MagneticMomentModule (Source104)
class MagneticMomentModuleTerm : public PhysicsTerm {
private:
    MagneticMomentModule instance;
public:
    MagneticMomentModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeMu_j(t);
    }
    
    std::string getName() const override {
        return "MagneticMomentModule";
    }
    
    std::string getDescription() const override {
        return "Source104: MagneticMomentModule.computeMu_j()";
    }
};

// Wrapper for GalacticBlackHoleModule (Source105)
class GalacticBlackHoleModuleTerm : public PhysicsTerm {
private:
    GalacticBlackHoleModule instance;
public:
    GalacticBlackHoleModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeM_bh(t);
    }
    
    std::string getName() const override {
        return "GalacticBlackHoleModule";
    }
    
    std::string getDescription() const override {
        return "Source105: GalacticBlackHoleModule.computeM_bh()";
    }
};

// Wrapper for NegativeTimeModule (Source106)
class NegativeTimeModuleTerm : public PhysicsTerm {
private:
    NegativeTimeModule instance;
public:
    NegativeTimeModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeT_n(t);
    }
    
    std::string getName() const override {
        return "NegativeTimeModule";
    }
    
    std::string getDescription() const override {
        return "Source106: NegativeTimeModule.computeT_n()";
    }
};

// Wrapper for PiConstantModule (Source107)
class PiConstantModuleTerm : public PhysicsTerm {
private:
    PiConstantModule instance;
public:
    PiConstantModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computePi(t);
    }
    
    std::string getName() const override {
        return "PiConstantModule";
    }
    
    std::string getDescription() const override {
        return "Source107: PiConstantModule.computePi()";
    }
};

// Wrapper for CorePenetrationModule (Source108)
class CorePenetrationModuleTerm : public PhysicsTerm {
private:
    CorePenetrationModule instance;
public:
    CorePenetrationModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeP_core(t);
    }
    
    std::string getName() const override {
        return "CorePenetrationModule";
    }
    
    std::string getDescription() const override {
        return "Source108: CorePenetrationModule.computeP_core()";
    }
};

// Wrapper for QuasiLongitudinalModule (Source109)
class QuasiLongitudinalModuleTerm : public PhysicsTerm {
private:
    QuasiLongitudinalModule instance;
public:
    QuasiLongitudinalModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF_quasi(t);
    }
    
    std::string getName() const override {
        return "QuasiLongitudinalModule";
    }
    
    std::string getDescription() const override {
        return "Source109: QuasiLongitudinalModule.computeF_quasi()";
    }
};

// Wrapper for OuterFieldBubbleModule (Source110)
class OuterFieldBubbleModuleTerm : public PhysicsTerm {
private:
    OuterFieldBubbleModule instance;
public:
    OuterFieldBubbleModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeR_b(t);
    }
    
    std::string getName() const override {
        return "OuterFieldBubbleModule";
    }
    
    std::string getDescription() const override {
        return "Source110: OuterFieldBubbleModule.computeR_b()";
    }
};

// Wrapper for ReciprocationDecayModule (Source111)
class ReciprocationDecayModuleTerm : public PhysicsTerm {
private:
    ReciprocationDecayModule instance;
public:
    ReciprocationDecayModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeGamma_day(t);
    }
    
    std::string getName() const override {
        return "ReciprocationDecayModule";
    }
    
    std::string getDescription() const override {
        return "Source111: ReciprocationDecayModule.computeGamma_day()";
    }
};

// Wrapper for ScmPenetrationModule (Source112)
class ScmPenetrationModuleTerm : public PhysicsTerm {
private:
    ScmPenetrationModule instance;
public:
    ScmPenetrationModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeP_SCm(t);
    }
    
    std::string getName() const override {
        return "ScmPenetrationModule";
    }
    
    std::string getDescription() const override {
        return "Source112: ScmPenetrationModule.computeP_SCm()";
    }
};

// Wrapper for ScmReactivityDecayModule (Source113)
class ScmReactivityDecayModuleTerm : public PhysicsTerm {
private:
    ScmReactivityDecayModule instance;
public:
    ScmReactivityDecayModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeKappa_day(t);
    }
    
    std::string getName() const override {
        return "ScmReactivityDecayModule";
    }
    
    std::string getDescription() const override {
        return "Source113: ScmReactivityDecayModule.computeKappa_day()";
    }
};

// Wrapper for SolarCycleFrequencyModule (Source114)
class SolarCycleFrequencyModuleTerm : public PhysicsTerm {
private:
    SolarCycleFrequencyModule instance;
public:
    SolarCycleFrequencyModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeOmega_c(t);
    }
    
    std::string getName() const override {
        return "SolarCycleFrequencyModule";
    }
    
    std::string getDescription() const override {
        return "Source114: SolarCycleFrequencyModule.computeOmega_c()";
    }
};

// Wrapper for SolarWindModulationModule (Source115)
class SolarWindModulationModuleTerm : public PhysicsTerm {
private:
    SolarWindModulationModule instance;
public:
    SolarWindModulationModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeDelta_sw(t);
    }
    
    std::string getName() const override {
        return "SolarWindModulationModule";
    }
    
    std::string getDescription() const override {
        return "Source115: SolarWindModulationModule.computeDelta_sw()";
    }
};

// Wrapper for SolarWindVelocityModule (Source116)
class SolarWindVelocityModuleTerm : public PhysicsTerm {
private:
    SolarWindVelocityModule instance;
public:
    SolarWindVelocityModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeV_sw(t);
    }
    
    std::string getName() const override {
        return "SolarWindVelocityModule";
    }
    
    std::string getDescription() const override {
        return "Source116: SolarWindVelocityModule.computeV_sw()";
    }
};

// Wrapper for StellarMassModule (Source117)
class StellarMassModuleTerm : public PhysicsTerm {
private:
    StellarMassModule instance;
public:
    StellarMassModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeM_s(t);
    }
    
    std::string getName() const override {
        return "StellarMassModule";
    }
    
    std::string getDescription() const override {
        return "Source117: StellarMassModule.computeM_s()";
    }
};

// Wrapper for StellarRotationModule (Source118)
class StellarRotationModuleTerm : public PhysicsTerm {
private:
    StellarRotationModule instance;
public:
    StellarRotationModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeOmega_s(t);
    }
    
    std::string getName() const override {
        return "StellarRotationModule";
    }
    
    std::string getDescription() const override {
        return "Source118: StellarRotationModule.computeOmega_s()";
    }
};

// Wrapper for StepFunctionModule (Source119)
class StepFunctionModuleTerm : public PhysicsTerm {
private:
    StepFunctionModule instance;
public:
    StepFunctionModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeS_r_Rb(t);
    }
    
    std::string getName() const override {
        return "StepFunctionModule";
    }
    
    std::string getDescription() const override {
        return "Source119: StepFunctionModule.computeS_r_Rb()";
    }
};

// Wrapper for StressEnergyTensorModule (Source120)
class StressEnergyTensorModuleTerm : public PhysicsTerm {
private:
    StressEnergyTensorModule instance;
public:
    StressEnergyTensorModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeT_s(t);
    }
    
    std::string getName() const override {
        return "StressEnergyTensorModule";
    }
    
    std::string getDescription() const override {
        return "Source120: StressEnergyTensorModule.computeT_s()";
    }
};

// Wrapper for SurfaceMagneticFieldModule (Source121)
class SurfaceMagneticFieldModuleTerm : public PhysicsTerm {
private:
    SurfaceMagneticFieldModule instance;
public:
    SurfaceMagneticFieldModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeB_s_min(t);
    }
    
    std::string getName() const override {
        return "SurfaceMagneticFieldModule";
    }
    
    std::string getDescription() const override {
        return "Source121: SurfaceMagneticFieldModule.computeB_s_min()";
    }
};

// Wrapper for SurfaceTemperatureModule (Source122)
class SurfaceTemperatureModuleTerm : public PhysicsTerm {
private:
    SurfaceTemperatureModule instance;
public:
    SurfaceTemperatureModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeT_s(t);
    }
    
    std::string getName() const override {
        return "SurfaceTemperatureModule";
    }
    
    std::string getDescription() const override {
        return "Source122: SurfaceTemperatureModule.computeT_s()";
    }
};

// Wrapper for TimeReversalZoneModule (Source123)
class TimeReversalZoneModuleTerm : public PhysicsTerm {
private:
    TimeReversalZoneModule instance;
public:
    TimeReversalZoneModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF_TRZ(t);
    }
    
    std::string getName() const override {
        return "TimeReversalZoneModule";
    }
    
    std::string getDescription() const override {
        return "Source123: TimeReversalZoneModule.computeF_TRZ()";
    }
};

// Wrapper for Ug1DefectModule (Source124)
class Ug1DefectModuleTerm : public PhysicsTerm {
private:
    Ug1DefectModule instance;
public:
    Ug1DefectModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeDelta_def(t);
    }
    
    std::string getName() const override {
        return "Ug1DefectModule";
    }
    
    std::string getDescription() const override {
        return "Source124: Ug1DefectModule.computeDelta_def()";
    }
};

// Wrapper for Ug3DiskVectorModule (Source125)
class Ug3DiskVectorModuleTerm : public PhysicsTerm {
private:
    Ug3DiskVectorModule instance;
public:
    Ug3DiskVectorModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computePhiHatMagnitude(t);
    }
    
    std::string getName() const override {
        return "Ug3DiskVectorModule";
    }
    
    std::string getDescription() const override {
        return "Source125: Ug3DiskVectorModule.computePhiHatMagnitude()";
    }
};

// Wrapper for AetherVacuumDensityModule (Source126)
class AetherVacuumDensityModuleTerm : public PhysicsTerm {
private:
    AetherVacuumDensityModule instance;
public:
    AetherVacuumDensityModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeRho_vac_A(t);
    }
    
    std::string getName() const override {
        return "AetherVacuumDensityModule";
    }
    
    std::string getDescription() const override {
        return "Source126: AetherVacuumDensityModule.computeRho_vac_A()";
    }
};

// Wrapper for UniversalInertiaVacuumModule (Source127)
class UniversalInertiaVacuumModuleTerm : public PhysicsTerm {
private:
    UniversalInertiaVacuumModule instance;
public:
    UniversalInertiaVacuumModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeRho_vac_Ui(t);
    }
    
    std::string getName() const override {
        return "UniversalInertiaVacuumModule";
    }
    
    std::string getDescription() const override {
        return "Source127: UniversalInertiaVacuumModule.computeRho_vac_Ui()";
    }
};

// Wrapper for ScmVacuumDensityModule (Source128)
class ScmVacuumDensityModuleTerm : public PhysicsTerm {
private:
    ScmVacuumDensityModule instance;
public:
    ScmVacuumDensityModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeRho_vac_SCm(t);
    }
    
    std::string getName() const override {
        return "ScmVacuumDensityModule";
    }
    
    std::string getDescription() const override {
        return "Source128: ScmVacuumDensityModule.computeRho_vac_SCm()";
    }
};

// Wrapper for UaVacuumDensityModule (Source129)
class UaVacuumDensityModuleTerm : public PhysicsTerm {
private:
    UaVacuumDensityModule instance;
public:
    UaVacuumDensityModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeRho_vac_UA(t);
    }
    
    std::string getName() const override {
        return "UaVacuumDensityModule";
    }
    
    std::string getDescription() const override {
        return "Source129: UaVacuumDensityModule.computeRho_vac_UA()";
    }
};

// Wrapper for UniversalInertiaVacuumModule (Source130)
class UniversalInertiaVacuumModuleTerm : public PhysicsTerm {
private:
    UniversalInertiaVacuumModule instance;
public:
    UniversalInertiaVacuumModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeRho_vac_Ui(t);
    }
    
    std::string getName() const override {
        return "UniversalInertiaVacuumModule";
    }
    
    std::string getDescription() const override {
        return "Source130: UniversalInertiaVacuumModule.computeRho_vac_Ui()";
    }
};

// Wrapper for ScmVelocityModule (Source131)
class ScmVelocityModuleTerm : public PhysicsTerm {
private:
    ScmVelocityModule instance;
public:
    ScmVelocityModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeE_react(t);
    }
    
    std::string getName() const override {
        return "ScmVelocityModule";
    }
    
    std::string getDescription() const override {
        return "Source131: ScmVelocityModule.computeE_react()";
    }
};

// Wrapper for ButterflyNebulaUQFFModule (Source132)
class ButterflyNebulaUQFFModuleTerm : public PhysicsTerm {
private:
    ButterflyNebulaUQFFModule instance;
public:
    ButterflyNebulaUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF_U_Bi(t);
    }
    
    std::string getName() const override {
        return "ButterflyNebulaUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source132: ButterflyNebulaUQFFModule.computeF_U_Bi()";
    }
};

// Wrapper for CentaurusAUQFFModule (Source133)
class CentaurusAUQFFModuleTerm : public PhysicsTerm {
private:
    CentaurusAUQFFModule instance;
public:
    CentaurusAUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF_U_Bi(t);
    }
    
    std::string getName() const override {
        return "CentaurusAUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source133: CentaurusAUQFFModule.computeF_U_Bi()";
    }
};

// Wrapper for Abell2256UQFFModule (Source134)
class Abell2256UQFFModuleTerm : public PhysicsTerm {
private:
    Abell2256UQFFModule instance;
public:
    Abell2256UQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF(t);
    }
    
    std::string getName() const override {
        return "Abell2256UQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source134: Abell2256UQFFModule.computeF()";
    }
};

// Wrapper for ASASSN14liUQFFModule (Source135)
class ASASSN14liUQFFModuleTerm : public PhysicsTerm {
private:
    ASASSN14liUQFFModule instance;
public:
    ASASSN14liUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF(t);
    }
    
    std::string getName() const override {
        return "ASASSN14liUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source135: ASASSN14liUQFFModule.computeF()";
    }
};

// Wrapper for CentaurusAUQFFModule (Source136)
class CentaurusAUQFFModuleTerm : public PhysicsTerm {
private:
    CentaurusAUQFFModule instance;
public:
    CentaurusAUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF(t);
    }
    
    std::string getName() const override {
        return "CentaurusAUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source136: CentaurusAUQFFModule.computeF()";
    }
};

// Wrapper for CrabNebulaUQFFModule (Source137)
class CrabNebulaUQFFModuleTerm : public PhysicsTerm {
private:
    CrabNebulaUQFFModule instance;
public:
    CrabNebulaUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF(t);
    }
    
    std::string getName() const override {
        return "CrabNebulaUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source137: CrabNebulaUQFFModule.computeF()";
    }
};

// Wrapper for ElGordoUQFFModule (Source138)
class ElGordoUQFFModuleTerm : public PhysicsTerm {
private:
    ElGordoUQFFModule instance;
public:
    ElGordoUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF(t);
    }
    
    std::string getName() const override {
        return "ElGordoUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source138: ElGordoUQFFModule.computeF()";
    }
};

// Wrapper for ESO137UQFFModule (Source139)
class ESO137UQFFModuleTerm : public PhysicsTerm {
private:
    ESO137UQFFModule instance;
public:
    ESO137UQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF(t);
    }
    
    std::string getName() const override {
        return "ESO137UQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source139: ESO137UQFFModule.computeF()";
    }
};

// Wrapper for IC2163UQFFModule (Source140)
class IC2163UQFFModuleTerm : public PhysicsTerm {
private:
    IC2163UQFFModule instance;
public:
    IC2163UQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF(t);
    }
    
    std::string getName() const override {
        return "IC2163UQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source140: IC2163UQFFModule.computeF()";
    }
};

// Wrapper for J1610UQFFModule (Source141)
class J1610UQFFModuleTerm : public PhysicsTerm {
private:
    J1610UQFFModule instance;
public:
    J1610UQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF(t);
    }
    
    std::string getName() const override {
        return "J1610UQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source141: J1610UQFFModule.computeF()";
    }
};

// Wrapper for JupiterAuroraeUQFFModule (Source142)
class JupiterAuroraeUQFFModuleTerm : public PhysicsTerm {
private:
    JupiterAuroraeUQFFModule instance;
public:
    JupiterAuroraeUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF(t);
    }
    
    std::string getName() const override {
        return "JupiterAuroraeUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source142: JupiterAuroraeUQFFModule.computeF()";
    }
};

// Wrapper for LagoonNebulaUQFFModule (Source143)
class LagoonNebulaUQFFModuleTerm : public PhysicsTerm {
private:
    LagoonNebulaUQFFModule instance;
public:
    LagoonNebulaUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF(t);
    }
    
    std::string getName() const override {
        return "LagoonNebulaUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source143: LagoonNebulaUQFFModule.computeF()";
    }
};

// Wrapper for LagoonNebulaUQFFModule (Source144)
class LagoonNebulaUQFFModuleTerm : public PhysicsTerm {
private:
    LagoonNebulaUQFFModule instance;
public:
    LagoonNebulaUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF(t);
    }
    
    std::string getName() const override {
        return "LagoonNebulaUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source144: LagoonNebulaUQFFModule.computeF()";
    }
};

// Wrapper for M87JetUQFFModule (Source145)
class M87JetUQFFModuleTerm : public PhysicsTerm {
private:
    M87JetUQFFModule instance;
public:
    M87JetUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF(t);
    }
    
    std::string getName() const override {
        return "M87JetUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source145: M87JetUQFFModule.computeF()";
    }
};

// Wrapper for NGC1365UQFFModule (Source146)
class NGC1365UQFFModuleTerm : public PhysicsTerm {
private:
    NGC1365UQFFModule instance;
public:
    NGC1365UQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF(t);
    }
    
    std::string getName() const override {
        return "NGC1365UQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source146: NGC1365UQFFModule.computeF()";
    }
};

// Wrapper for NGC2207UQFFModule (Source147)
class NGC2207UQFFModuleTerm : public PhysicsTerm {
private:
    NGC2207UQFFModule instance;
public:
    NGC2207UQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF(t);
    }
    
    std::string getName() const override {
        return "NGC2207UQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source147: NGC2207UQFFModule.computeF()";
    }
};

// Wrapper for RAquariiUQFFModule (Source148)
class RAquariiUQFFModuleTerm : public PhysicsTerm {
private:
    RAquariiUQFFModule instance;
public:
    RAquariiUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF(t);
    }
    
    std::string getName() const override {
        return "RAquariiUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source148: RAquariiUQFFModule.computeF()";
    }
};

// Wrapper for SgrAStarUQFFModule (Source149)
class SgrAStarUQFFModuleTerm : public PhysicsTerm {
private:
    SgrAStarUQFFModule instance;
public:
    SgrAStarUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF(t);
    }
    
    std::string getName() const override {
        return "SgrAStarUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source149: SgrAStarUQFFModule.computeF()";
    }
};

// Wrapper for SPTCLJ2215UQFFModule (Source150)
class SPTCLJ2215UQFFModuleTerm : public PhysicsTerm {
private:
    SPTCLJ2215UQFFModule instance;
public:
    SPTCLJ2215UQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF(t);
    }
    
    std::string getName() const override {
        return "SPTCLJ2215UQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source150: SPTCLJ2215UQFFModule.computeF()";
    }
};

// Wrapper for StephanQuintetUQFFModule (Source151)
class StephanQuintetUQFFModuleTerm : public PhysicsTerm {
private:
    StephanQuintetUQFFModule instance;
public:
    StephanQuintetUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF(t);
    }
    
    std::string getName() const override {
        return "StephanQuintetUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source151: StephanQuintetUQFFModule.computeF()";
    }
};

// Wrapper for VelaPulsarUQFFModule (Source152)
class VelaPulsarUQFFModuleTerm : public PhysicsTerm {
private:
    VelaPulsarUQFFModule instance;
public:
    VelaPulsarUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF(t);
    }
    
    std::string getName() const override {
        return "VelaPulsarUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source152: VelaPulsarUQFFModule.computeF()";
    }
};

// Wrapper for Abell2256UQFFModule (Source153)
class Abell2256UQFFModuleTerm : public PhysicsTerm {
private:
    Abell2256UQFFModule instance;
public:
    Abell2256UQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeF(t);
    }
    
    std::string getName() const override {
        return "Abell2256UQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source153: Abell2256UQFFModule.computeF()";
    }
};

// Wrapper for HydrogenResonanceUQFFModule (Source154)
class HydrogenResonanceUQFFModuleTerm : public PhysicsTerm {
private:
    HydrogenResonanceUQFFModule instance;
public:
    HydrogenResonanceUQFFModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeHRes(t);
    }
    
    std::string getName() const override {
        return "HydrogenResonanceUQFFModule";
    }
    
    std::string getDescription() const override {
        return "Source154: HydrogenResonanceUQFFModule.computeHRes()";
    }
};

// Wrapper for UQFFBuoyancyModule (Source155)
class UQFFBuoyancyModuleTerm : public PhysicsTerm {
private:
    UQFFBuoyancyModule instance;
public:
    UQFFBuoyancyModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeFBi(t);
    }
    
    std::string getName() const override {
        return "UQFFBuoyancyModule";
    }
    
    std::string getDescription() const override {
        return "Source155: UQFFBuoyancyModule.computeFBi()";
    }
};

// Wrapper for UQFFBuoyancyCNBModule (Source156)
class UQFFBuoyancyCNBModuleTerm : public PhysicsTerm {
private:
    UQFFBuoyancyCNBModule instance;
public:
    UQFFBuoyancyCNBModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeFBi(t);
    }
    
    std::string getName() const override {
        return "UQFFBuoyancyCNBModule";
    }
    
    std::string getDescription() const override {
        return "Source156: UQFFBuoyancyCNBModule.computeFBi()";
    }
};

// Wrapper for SurfaceMagneticFieldModule (Source157)
class SurfaceMagneticFieldModuleTerm : public PhysicsTerm {
private:
    SurfaceMagneticFieldModule instance;
public:
    SurfaceMagneticFieldModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeUb1(t);
    }
    
    std::string getName() const override {
        return "SurfaceMagneticFieldModule";
    }
    
    std::string getDescription() const override {
        return "Source157: SurfaceMagneticFieldModule.computeUb1()";
    }
};

// Wrapper for UQFFBuoyancyModule (Source158)
class UQFFBuoyancyModuleTerm : public PhysicsTerm {
private:
    UQFFBuoyancyModule instance;
public:
    UQFFBuoyancyModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeFBi(t);
    }
    
    std::string getName() const override {
        return "UQFFBuoyancyModule";
    }
    
    std::string getDescription() const override {
        return "Source158: UQFFBuoyancyModule.computeFBi()";
    }
};

// Wrapper for UQFFBuoyancyModule (Source159)
class UQFFBuoyancyModuleTerm : public PhysicsTerm {
private:
    UQFFBuoyancyModule instance;
public:
    UQFFBuoyancyModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeFBi(t);
    }
    
    std::string getName() const override {
        return "UQFFBuoyancyModule";
    }
    
    std::string getDescription() const override {
        return "Source159: UQFFBuoyancyModule.computeFBi()";
    }
};

// Wrapper for UQFFBuoyancyModule (Source160)
class UQFFBuoyancyModuleTerm : public PhysicsTerm {
private:
    UQFFBuoyancyModule instance;
public:
    UQFFBuoyancyModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeFBi(t);
    }
    
    std::string getName() const override {
        return "UQFFBuoyancyModule";
    }
    
    std::string getDescription() const override {
        return "Source160: UQFFBuoyancyModule.computeFBi()";
    }
};

// Wrapper for UQFFBuoyancyAstroModule (Source161)
class UQFFBuoyancyAstroModuleTerm : public PhysicsTerm {
private:
    UQFFBuoyancyAstroModule instance;
public:
    UQFFBuoyancyAstroModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeFBi(t);
    }
    
    std::string getName() const override {
        return "UQFFBuoyancyAstroModule";
    }
    
    std::string getDescription() const override {
        return "Source161: UQFFBuoyancyAstroModule.computeFBi()";
    }
};

// Wrapper for UQFFBuoyancyCNBModule (Source162)
class UQFFBuoyancyCNBModuleTerm : public PhysicsTerm {
private:
    UQFFBuoyancyCNBModule instance;
public:
    UQFFBuoyancyCNBModuleTerm() {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Call embedded module's actual compute method
        return instance.computeFBi(t);
    }
    
    std::string getName() const override {
        return "UQFFBuoyancyCNBModule";
    }
    
    std::string getDescription() const override {
        return "Source162: UQFFBuoyancyCNBModule.computeFBi()";
    }
};

