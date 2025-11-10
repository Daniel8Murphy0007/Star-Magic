# CoAnQi Quick Reference

## Compilation & Execution

```bash
g++ -std=c++17 MAIN_1_CoAnQi.cpp -o MAIN_1_CoAnQi
./MAIN_1_CoAnQi
```

## Main Menu Quick Reference

| Option | Function | What It Does |
|--------|----------|--------------|
| 1 | Calculate System | Compute UQFF for single system |
| 2 | Calculate ALL Systems | Batch process with statistics |
| 3 | Clone & Mutate | Generate derivative system |
| 4 | Add Custom System | Define new system |
| 5 | Add Physics Term | Generate term code |
| 6 | Run Simulations | Execute 6 simulation modes |
| 7 | Statistical Analysis | Full statistical suite |
| 8 | Self-Optimization | Auto-tune parameters |
| 9 | Exit | Shutdown |

## Key Capabilities

✅ **Self-Expanding**: Add physics terms at runtime  
✅ **Self-Updating**: Optimizes parameters automatically  
✅ **Self-Cloning**: Generates derivative systems  
✅ **Statistical**: Full analysis suite (mean, stddev, correlation)  
✅ **Verbose**: 3-level logging (INFO/CALC/DEBUG)  
✅ **Simultaneous**: Process all systems at once  

## Physics Terms Integrated

1. **DynamicVacuumTerm** - Time-varying vacuum energy
2. **QuantumCouplingTerm** - Quantum entanglement
3. **DarkMatterHaloTerm** - NFW halo contributions
4. **VacuumEnergyTerm** - Large-scale vacuum variations
5. **QuantumEntanglementTerm** - Spooky action effects
6. **CosmicNeutrinoTerm** - CNB contributions

## Core UQFF Outputs

- **F_U_Bi_i**: UQFF buoyancy force (N)
- **g_compressed**: 26-layer gravity field (m/s²)
- **Dynamic terms**: Additional physics (N)
- **F_jet_rel**: Relativistic jet thrust (N)
- **E_acc_rel**: Acceleration coherence (J)
- **F_drag_rel**: Relativistic drag (N)
- **F_gw_rel**: GW ripple force (N)

## Example Workflow

```
1. Start CoAnQi
2. Choose option 1 (single system)
3. Enter "Vela Pulsar"
4. View comprehensive results
5. Choose option 3 (clone)
6. Enter "Vela Pulsar", mutation 0.1
7. Choose option 2 (all systems)
8. View statistical analysis
9. Exit
```

## Log Files

- **Location**: `coAnQi_log_<timestamp>.txt`
- **Levels**: 1=INFO, 2=CALC, 3=DEBUG
- **Format**: `<timestamp> [LEVEL] <message>`

## Comparison: MAIN_1.cpp vs MAIN_1_CoAnQi.cpp

| Feature | MAIN_1.cpp | MAIN_1_CoAnQi.cpp |
|---------|-----------|-------------------|
| Lines of Code | 1720 | ~1500 |
| Physics Terms | Core UQFF only | Core + 6 dynamic terms |
| Self-Expanding | ❌ No | ✅ Yes |
| Self-Updating | ❌ No | ✅ Yes |
| Self-Cloning | ❌ No | ✅ Yes |
| Statistical Analysis | ❌ No | ✅ Yes |
| Verbose Logging | ❌ No | ✅ Yes |
| Batch Processing | ❌ No | ✅ Yes |
| Code Generation | ❌ No | ✅ Yes |
| Module Registry | ❌ No | ✅ Yes |
| Parameter Optimization | ❌ No | ✅ Yes |

## Next Steps

1. **Explore**: Run single systems, observe outputs
2. **Clone**: Create derivative systems with mutations
3. **Analyze**: Use statistical tools on batch results
4. **Optimize**: Let system auto-tune parameters
5. **Extend**: Add custom physics terms
6. **Integrate**: Use alongside Source13-162 modules

---
**CoAnQi = Conscious Quantum Intelligence**  
*Where physics becomes alive!* 🌟
