# UQFF Star-Magic — Commercial License Information

Copyright (c) 2025-2026 Daniel T. Murphy / Star-Magic Research Program.
All rights reserved.

This document explains when you need a commercial license for the UQFF
Star-Magic codebase, what the commercial license typically covers, and
how to request one. The default open-source option is AGPL-3.0
(see `LICENSE` and `LICENSE-AGPL-3.0.txt`).

---

## 1. When do I need a commercial license?

You can use UQFF Star-Magic for free under AGPL-3.0 for:

- Academic teaching and research
- Replication studies and peer review
- Personal exploration
- Non-commercial open-source projects that themselves accept AGPL-3.0

You need a **commercial license** if any of these apply:

| Scenario | AGPL ok? | Commercial license required? |
|---|---|---|
| University researcher publishes a paper using UQFF | Yes | No |
| Private lab uses UQFF internally for R&D, no distribution | Yes | No |
| Software vendor ships proprietary product that links UQFF | No | **Yes** |
| Cloud company offers UQFF-as-a-Service to paying users | No | **Yes** |
| Reactor manufacturer embeds UQFF in firmware sold commercially | No | **Yes** |
| Energy utility uses UQFF predictions to operate a commercial plant | Depends | **Likely yes** |
| Government national lab uses UQFF for classified work | Negotiable | **Probably yes** |
| University spin-off commercializes UQFF-derived patents | No | **Yes** |
| Open-source MIT-licensed project wants to import UQFF | No | **Yes** (or relicense your project to AGPL) |

The key triggers are:
1. You want to keep your derivative source **closed**.
2. You want to **distribute or host** UQFF without inheriting AGPL.
3. You're commercializing **Star-Magic hardware** (LENR reactors, etc.).

---

## 2. What does a commercial license typically cover?

A standard UQFF Star-Magic commercial license grants:

- A perpetual or term-limited right to use, modify, and distribute the
  software in your proprietary products without AGPL obligations.
- A patent license covering relevant Star-Magic patents for the
  agreed field of use.
- Optional source-modification rights (for keeping internal forks).
- Optional support and consultation hours with the author.
- Optional indemnification against third-party IP claims.
- Optional field-of-use exclusivity (subject to upfront premium).

Terms are negotiated case-by-case based on:

- Your **field of use** (research / SaaS / hardware / consulting / etc.)
- Your **deployment scale** (number of users, number of installations,
  number of CPU cores, revenue from product, etc.)
- Whether you need **source-modification rights**.
- Whether you need **support and indemnification**.
- **Duration** (perpetual, annual, multi-year).
- **Exclusivity** (non-exclusive is default; exclusive is expensive).

---

## 3. How do I request a commercial license?

Email: **daniel.murphy00@enrgyone.com**
Subject: **"UQFF Star-Magic Commercial License Request"**

Please include the following:

```
Organization name:
Country of incorporation:
Primary contact name + title:
Intended field of use:
  [ ] Internal research only
  [ ] Distributed software product
  [ ] SaaS / web-API
  [ ] Embedded firmware (reactor / sensor / instrument)
  [ ] Consulting / advisory
  [ ] Other: ____________
Expected deployment scale (users, installs, revenue):
Duration desired (1 year / 3 years / perpetual):
Source-modification rights needed:  [ ] Yes  [ ] No
Support hours desired (0 / 10 / 50 / 100+):
Indemnification needed:             [ ] Yes  [ ] No
Field-of-use exclusivity desired:   [ ] Non-exclusive  [ ] Exclusive
Target start date:
Brief description of your product / use case:
```

You will receive a response within 10 business days with either:
- A draft term sheet for negotiation, or
- A request for additional information, or
- A polite decline (e.g., conflict with another exclusive licensee).

---

## 4. Frequently Asked Questions

### "Can I evaluate UQFF Star-Magic before paying?"

Yes. The AGPL-3.0 option lets you evaluate the full codebase for any
length of time. The trigger for needing a commercial license is
**distribution or hosted-service deployment**, not evaluation.

### "What if I'm a small startup?"

There is a startup-friendly tier for companies with <$5M annual
revenue, <20 employees, and <3 years since incorporation. Mention
your stage in your request and we'll discuss reduced terms.

### "What if I'm an individual contractor?"

If you're an individual using UQFF on behalf of a client, the client
needs the license, not you. Refer the client to this document.

### "Does my open-source project need a commercial license?"

If your project is itself AGPL-3.0 (or GPL-3.0-compatible), no — you
can use UQFF freely. If your project is MIT, Apache-2.0, BSD, or any
other license that's incompatible with AGPL-3.0's share-alike clause,
you either need a commercial license OR you need to relicense your
project to AGPL-3.0.

### "Can I publish papers that cite UQFF without a commercial license?"

Yes. Academic publication is explicitly permitted under AGPL-3.0 and
requires no commercial license. You MUST cite UQFF per `CITATION.cff`.

### "Can my university patent-protect work that uses UQFF?"

Yes for the patent itself. However, if you want to **commercialize**
the patent (e.g., license it to industry), the commercial product
incorporating UQFF will need a commercial license. Talk to your
technology transfer office and contact us early.

### "What if I need to sublicense UQFF to my customers?"

Sublicensing rights are available but require negotiation. The default
commercial license is single-organization; sublicense rights add a
premium and require revenue reporting.

### "What about military / defense / classified use?"

Defense and classified use is negotiable but requires additional
review including (in some jurisdictions) export-control clearance.
Contact us before assuming AGPL-3.0 alone is sufficient.

### "What jurisdiction governs the commercial license?"

The default governing law is the jurisdiction of incorporation of the
licensor (currently the United States, individual sole-trader). This
is negotiable for large licensees.

---

## 5. Trademark and Patent Notes

- "UQFF", "Star-Magic", "Di-Pseudo-Monopole", and "DPM" are
  unregistered trademarks of Daniel T. Murphy / Star-Magic Research
  Program. Neither AGPL-3.0 nor a commercial license grants trademark
  rights. You may factually state your product "uses UQFF" but may
  not use the marks in any manner implying endorsement.

- The Star-Magic LENR reactor architecture (COP 555:1 at 27 W,
  ambient T, pH -37) may be subject to pending or issued patents.
  Software licenses do NOT grant hardware patent rights. A separate
  patent license is required for hardware implementations.

---

End of commercial license information.
