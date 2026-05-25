#!/usr/bin/env python3
"""
ENTANGLEMENT TRANSACTION LEDGER - Tier 1 Foundation

Pillar 2: Detailed tracking of energy lending between electron twins

Physical Process:
- Twin electrons share entanglement bond
- Electron 1 can "borrow" energy from electron 2 (and vice versa)
- Lending tracks accumulated imbalance
- Repayment happens when entanglement oscillates

Key Insight: This explains how atoms maintain coherence despite energy fluctuations
- No external energy needed
- Atoms self-regulate through internal energy redistribution
- Lending capacity ~DPM binding energy

Mathematical Reference: COMPLETE_UQFF_UNIFIED_FRAMEWORK.md (Part II.3)
Date: May 24, 2026
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
from collections import deque
from enum import Enum

# ============================================================================
# ENUM: TRANSACTION TYPES
# ============================================================================

class TransactionType(Enum):
    """Types of energy transactions between twins"""
    BORROW = "borrow"          # Electron 1 borrows from electron 2
    REPAY = "repay"            # Electron 1 repays to electron 2
    OSCILLATE = "oscillate"    # DPM field returns energy
    DECAY = "decay"            # Coherence time depletes lending capacity


# ============================================================================
# PHYSICAL CONSTANTS
# ============================================================================

@dataclass
class LedgerConstants:
    """Constants for transaction tracking"""
    HBAR: float = 1.054571817e-34
    M_E: float = 9.1093837e-31
    C: float = 299792458.0
    ALPHA_FINE: float = 1/137.036


# ============================================================================
# TRANSACTION RECORD
# ============================================================================

@dataclass
class EnergyTransaction:
    """Single energy transaction between twin electrons"""
    timestamp: float              # Time (s)
    transaction_type: TransactionType
    amount_ev: float             # Energy amount (eV)
    from_electron: int           # Source (1 or 2)
    to_electron: int             # Destination (1 or 2)
    reason: str                  # Physical explanation
    coherence_factor: float      # DPM coherence at transaction time (0-1)
    
    def __str__(self) -> str:
        return (f"T={self.timestamp:.3e}s | {self.transaction_type.value.upper()} "
                f"{self.amount_ev:.3e} eV | "
                f"e{self.from_electron}→e{self.to_electron} "
                f"(coherence={self.coherence_factor:.3f})")


# ============================================================================
# ELECTRON TWIN LEDGER
# ============================================================================

@dataclass
class TwinElectronLedger:
    """
    Complete transaction ledger for electron twin pair
    """
    Z: int                                           # Atomic number
    n: int                                           # Shell number
    separation_m: float                              # Spatial separation (m)
    
    # Balances
    electron_1_balance_ev: float = 0.0              # Net energy of electron 1
    electron_2_balance_ev: float = 0.0              # Net energy of electron 2
    
    # Transaction history
    transactions: List[EnergyTransaction] = field(default_factory=list)
    
    # Maximum lending capacity (set at initialization)
    max_lending_capacity_ev: float = 1.022e6        # Default: pair creation cost
    
    # State tracking
    total_borrowed_lifetime_ev: float = 0.0         # Cumulative borrowed
    total_repaid_lifetime_ev: float = 0.0           # Cumulative repaid
    lending_efficiency: float = 1.0                 # Repayment rate / borrow rate
    
    const: LedgerConstants = field(default_factory=LedgerConstants)
    
    # Rolling window for recent activity
    recent_transactions: deque = field(default_factory=lambda: deque(maxlen=100))
    
    def current_imbalance_ev(self) -> float:
        """
        Current energy imbalance between twins.
        
        Returns: |E1 - E2| (should stay << max_lending_capacity)
        """
        return abs(self.electron_1_balance_ev - self.electron_2_balance_ev)
    
    def available_to_borrow_ev(self) -> float:
        """
        Maximum energy electron 1 can currently borrow from electron 2.
        
        Formula: max_available = max_lending_capacity * coherence_strength(t)
        
        Returns: Energy in eV
        """
        # For now, assume high coherence (will be updated by caller with time-dependent coherence)
        return self.max_lending_capacity_ev * 0.9
    
    def borrow(self, amount_ev: float, timestamp: float, coherence_factor: float,
               reason: str = "orbital adjustment") -> Tuple[bool, str]:
        """
        Electron 1 borrows energy from electron 2.
        
        Args:
            amount_ev: Amount to borrow (eV)
            timestamp: Current time (s)
            coherence_factor: DPM coherence strength (0-1)
            reason: Why borrowing is needed
        
        Returns:
            (success, message)
        """
        available = self.max_lending_capacity_ev * coherence_factor
        
        if amount_ev > available:
            return (False, f"Insufficient lending capacity: "
                          f"requested {amount_ev:.3e} eV, "
                          f"available {available:.3e} eV")
        
        # Record transaction
        transaction = EnergyTransaction(
            timestamp=timestamp,
            transaction_type=TransactionType.BORROW,
            amount_ev=amount_ev,
            from_electron=2,
            to_electron=1,
            reason=reason,
            coherence_factor=coherence_factor
        )
        
        # Update balances
        self.electron_1_balance_ev += amount_ev
        self.electron_2_balance_ev -= amount_ev
        
        # Update lifetime counters
        self.total_borrowed_lifetime_ev += amount_ev
        
        # Record
        self.transactions.append(transaction)
        self.recent_transactions.append(transaction)
        
        return (True, f"Borrowed {amount_ev:.3e} eV (new balance: e1={self.electron_1_balance_ev:.3e}, "
                     f"e2={self.electron_2_balance_ev:.3e})")
    
    def repay(self, amount_ev: float, timestamp: float, reason: str = "DPM oscillation") -> Tuple[bool, str]:
        """
        Electron 1 repays borrowed energy back to electron 2.
        
        Args:
            amount_ev: Amount to repay (eV)
            timestamp: Current time (s)
            reason: Why repayment is happening
        
        Returns:
            (success, message)
        """
        if self.electron_1_balance_ev < amount_ev:
            return (False, f"Cannot repay: electron 1 balance {self.electron_1_balance_ev:.3e} eV "
                          f"< repayment {amount_ev:.3e} eV")
        
        # Record transaction
        transaction = EnergyTransaction(
            timestamp=timestamp,
            transaction_type=TransactionType.REPAY,
            amount_ev=amount_ev,
            from_electron=1,
            to_electron=2,
            reason=reason,
            coherence_factor=1.0  # Repayment assumed at good coherence
        )
        
        # Update balances
        self.electron_1_balance_ev -= amount_ev
        self.electron_2_balance_ev += amount_ev
        
        # Update lifetime counters
        self.total_repaid_lifetime_ev += amount_ev
        
        # Record
        self.transactions.append(transaction)
        self.recent_transactions.append(transaction)
        
        # Update lending efficiency
        if self.total_borrowed_lifetime_ev > 0:
            self.lending_efficiency = self.total_repaid_lifetime_ev / self.total_borrowed_lifetime_ev
        
        return (True, f"Repaid {amount_ev:.3e} eV (new balance: e1={self.electron_1_balance_ev:.3e}, "
                     f"e2={self.electron_2_balance_ev:.3e})")
    
    def sustainability_analysis(self) -> Dict:
        """
        Analyze lending sustainability.
        
        Checks:
        1. Is imbalance growing (signs of instability)?
        2. Is repayment keeping up with borrowing?
        3. Do transaction patterns suggest oscillation?
        
        Returns:
            Dictionary with sustainability metrics
        """
        imbalance = self.current_imbalance_ev()
        
        # Recent trend analysis
        n_recent = len(self.recent_transactions)
        if n_recent > 10:
            recent_subset = list(self.recent_transactions)[-10:]
            borrow_count = sum(1 for t in recent_subset if t.transaction_type == TransactionType.BORROW)
            repay_count = sum(1 for t in recent_subset if t.transaction_type == TransactionType.REPAY)
        else:
            borrow_count = repay_count = 0
        
        # Oscillation check: should see alternating borrow/repay
        has_oscillation = (borrow_count > 0) and (repay_count > 0)
        
        # Sustainability: imbalance should stay << max_capacity
        is_sustainable = imbalance < (0.5 * self.max_lending_capacity_ev)
        
        return {
            'current_imbalance_ev': imbalance,
            'total_borrowed_ev': self.total_borrowed_lifetime_ev,
            'total_repaid_ev': self.total_repaid_lifetime_ev,
            'lending_efficiency': self.lending_efficiency,
            'recent_borrow_count': borrow_count,
            'recent_repay_count': repay_count,
            'has_oscillation': has_oscillation,
            'is_sustainable': is_sustainable,
            'sustainability_status': "✓ STABLE" if is_sustainable else "✗ UNSTABLE",
        }
    
    def transaction_summary(self) -> str:
        """
        Summary of all transactions.
        
        Returns:
            Formatted summary string
        """
        lines = [
            f"LEDGER: Z={self.Z}, n={self.n}, separation={self.separation_m:.3e} m",
            f"  Total transactions: {len(self.transactions)}",
            f"  Electron 1 balance: {self.electron_1_balance_ev:.3e} eV",
            f"  Electron 2 balance: {self.electron_2_balance_ev:.3e} eV",
            f"  Imbalance: {self.current_imbalance_ev():.3e} eV",
            f"  Total borrowed (lifetime): {self.total_borrowed_lifetime_ev:.3e} eV",
            f"  Total repaid (lifetime): {self.total_repaid_lifetime_ev:.3e} eV",
            f"  Lending efficiency: {self.lending_efficiency:.3f}",
            f"  Max lending capacity: {self.max_lending_capacity_ev:.3e} eV",
        ]
        
        return "\n".join(lines)
    
    def transaction_log(self, n_recent: int = 20) -> str:
        """
        Display recent transaction log.
        
        Args:
            n_recent: Number of recent transactions to show
        
        Returns:
            Formatted transaction log
        """
        if not self.transactions:
            return "No transactions recorded."
        
        recent = self.transactions[-n_recent:] if n_recent else self.transactions
        
        lines = [f"TRANSACTION LOG (last {min(n_recent, len(recent))} of {len(self.transactions)}):"]
        for t in recent:
            lines.append(f"  {t}")
        
        return "\n".join(lines)


# ============================================================================
# LEDGER MANAGER
# ============================================================================

class LedgerManager:
    """
    Manage multiple twin pair ledgers across atomic system
    """
    
    def __init__(self):
        self.ledgers: Dict[Tuple[int, int], TwinElectronLedger] = {}
    
    def create_ledger(self, Z: int, n: int, separation_m: float,
                     max_lending_ev: float = 1.022e6) -> TwinElectronLedger:
        """
        Create new ledger for electron pair.
        """
        key = (Z, n)
        ledger = TwinElectronLedger(
            Z=Z,
            n=n,
            separation_m=separation_m,
            max_lending_capacity_ev=max_lending_ev
        )
        self.ledgers[key] = ledger
        return ledger
    
    def get_ledger(self, Z: int, n: int) -> Optional[TwinElectronLedger]:
        """Get existing ledger."""
        return self.ledgers.get((Z, n))
    
    def all_sustainable(self) -> bool:
        """Check if all ledgers are sustainable."""
        return all(ledger.sustainability_analysis()['is_sustainable'] 
                  for ledger in self.ledgers.values())
    
    def global_summary(self) -> str:
        """Summary across all ledgers."""
        lines = ["GLOBAL LEDGER SUMMARY", "=" * 60]
        
        total_transactions = sum(len(l.transactions) for l in self.ledgers.values())
        lines.append(f"Total ledgers: {len(self.ledgers)}")
        lines.append(f"Total transactions: {total_transactions}")
        
        total_borrowed = sum(l.total_borrowed_lifetime_ev for l in self.ledgers.values())
        total_repaid = sum(l.total_repaid_lifetime_ev for l in self.ledgers.values())
        lines.append(f"Total borrowed (system): {total_borrowed:.3e} eV")
        lines.append(f"Total repaid (system): {total_repaid:.3e} eV")
        
        is_sustainable = self.all_sustainable()
        lines.append(f"\nSystem sustainability: {'✓ STABLE' if is_sustainable else '✗ UNSTABLE'}")
        
        return "\n".join(lines)


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    print("=" * 80)
    print("ENTANGLEMENT TRANSACTION LEDGER - Tier 1 Test")
    print("=" * 80)
    
    # Create manager and test ledgers
    manager = LedgerManager()
    
    # Test case: Helium atom (2 electrons)
    ledger = manager.create_ledger(Z=2, n=1, separation_m=2e-12, 
                                  max_lending_ev=1.022e6)
    
    print("\nINITIAL STATE:")
    print(ledger.transaction_summary())
    
    # Simulate energy exchanges
    print("\n\nSIMULATED ENERGY EXCHANGES:\n")
    
    # Transaction 1: Borrow
    success, msg = ledger.borrow(1000, timestamp=0.0, coherence_factor=0.95,
                                reason="Orbital perturbation")
    print(f"1. Borrow: {msg}")
    
    # Transaction 2: Partial repay
    success, msg = ledger.repay(600, timestamp=1e-16, reason="DPM oscillation cycle 1")
    print(f"2. Repay: {msg}")
    
    # Transaction 3: Borrow again
    success, msg = ledger.borrow(400, timestamp=2e-16, coherence_factor=0.92,
                                reason="Fine structure oscillation")
    print(f"3. Borrow: {msg}")
    
    # Transaction 4: Full repay
    success, msg = ledger.repay(ledger.electron_1_balance_ev, timestamp=3e-16,
                               reason="DPM oscillation cycle 2")
    print(f"4. Repay all: {msg}")
    
    # Analysis
    print("\n\nSUSTAINABILITY ANALYSIS:")
    analysis = ledger.sustainability_analysis()
    for key, value in analysis.items():
        print(f"  {key}: {value}")
    
    # Final summary
    print("\n\nFINAL LEDGER STATE:")
    print(ledger.transaction_summary())
    
    print("\n\nTRANSACTION LOG:")
    print(ledger.transaction_log())
    
    print("\n" + "=" * 80)
    print("NEXT: Ledgers feed into simultaneous_7layer_solver.cpp")
    print("=" * 80)
