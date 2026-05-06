// makeMeasuredCartesianMap.C
//
// Writes the CERN-measured sPHENIX solenoid field map as a ROOT TNtuple
// on the same 111³ Cartesian grid used by PHField3DCartesian in offline
// tracking.  The output is a drop-in replacement for the OPERA map.
//
// Output TNtuple "fieldmap" branches: x, y, z (cm), bx, by, bz, hz (T)
// Grid: x/y/z in [-110, +110] cm, 2 cm steps, x slowest / z fastest —
//   identical layout to the original OPERA file.
//
// Advantages over the OPERA map:
//   - Ground-truth hardware measurement from CERN acceptance testing
//   - Satisfies ∇·B = 0 by construction (EnforceMaxwell in sPHENIXFieldMap)
//   - Correct axial position and field amplitude
//   - No spurious azimuthal dependence
//
// Points with r > 900 mm (outside measured coverage) are set to zero;
// these lie in the corners of the cube well outside the tracking volume.
//
// Run from the project directory:
//   root -l -b -q 'makeMeasuredCartesianMap.C+'
//
// Output: output/sphenix_measured_fieldmap_cartesian.root

#include "sphenix-cernfinal-map/sPHENIXFieldMap.cxx"

#include "TFile.h"
#include "TNtuple.h"
#include "TSystem.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

// ─────────────────────────────────────────────────────────────────────────────
void makeMeasuredCartesianMap()
{
    // ── Paths ─────────────────────────────────────────────────────────────────
    const char *FINE_CSV  =
        "/Users/haggerty/sphenix/data/data02/sphenix/MagnetMapping/cernfinal"
        "/fieldMapFineFullField.csv";
    const char *ROUGH_CSV =
        "/Users/haggerty/sphenix/data/data02/sphenix/MagnetMapping/cernfinal"
        "/fieldMapRoughFullField.csv";
    const char *OUT_FILE  =
        "output/sphenix_measured_fieldmap_cartesian.root";

    // ── Output grid (identical to OPERA analysis map) ─────────────────────────
    const int   N    = 111;       // nodes per axis
    const float GMIN = -110.f;    // cm
    const float DG   =    2.f;    // cm/step
    const float CM2MM = 10.f;     // unit conversion

    // ── Load measured field map ───────────────────────────────────────────────
    printf("=== Loading measured field map ===\n");
    sPHENIXFieldMap meas(FINE_CSV, ROUGH_CSV);
    printf("  Coverage: r = [%.0f, %.0f] mm,  z = [%.0f, %.0f] mm\n",
           meas.GetRMin(), meas.GetRMax(), meas.GetZMin(), meas.GetZMax());

    // ── Create output ─────────────────────────────────────────────────────────
    gSystem->mkdir("output", kTRUE);
    TFile *fOut = TFile::Open(OUT_FILE, "RECREATE");
    if (!fOut || fOut->IsZombie()) {
        fprintf(stderr, "ERROR: cannot create %s\n", OUT_FILE);
        return;
    }
    // Branch layout matches the original OPERA TNtuple exactly (hz = bz duplicate)
    TNtuple *nt = new TNtuple("fieldmap",
                              "sPHENIX solenoid field map — measured (Cartesian, cm/T)",
                              "x:y:z:bx:by:bz:hz");

    // ── Fill loop ─────────────────────────────────────────────────────────────
    // Traverse x slowest, y middle, z fastest — same order as original OPERA file.
    printf("=== Filling %d^3 = %d grid points ===\n", N, N*N*N);

    long nZero = 0;
    float vals[7];

    for (int ix = 0; ix < N; ++ix) {
        float x = GMIN + ix * DG;
        for (int iy = 0; iy < N; ++iy) {
            float y = GMIN + iy * DG;
            for (int iz = 0; iz < N; ++iz) {
                float z = GMIN + iz * DG;

                // GetFieldXYZ: position in mm, field in T
                double Bx, By, Bz;
                meas.GetFieldXYZ(x*CM2MM, y*CM2MM, z*CM2MM, Bx, By, Bz);

                vals[0] = x;
                vals[1] = y;
                vals[2] = z;
                vals[3] = (float)Bx;
                vals[4] = (float)By;
                vals[5] = (float)Bz;
                vals[6] = (float)Bz;   // hz = bz, matching original format
                nt->Fill(vals);

                if (Bz == 0. && Bx == 0. && By == 0.) ++nZero;
            }
        }
        if (ix % 11 == 0)
            printf("  ix = %3d / %d  (%.0f%%)\n", ix, N-1, 100.f*ix/(N-1));
    }
    printf("  %ld zero-field points (r > 900 mm or z outside measured range)\n",
           nZero);

    fOut->cd();
    nt->Write();
    fOut->Close();
    printf("  -> %s  (%d entries)\n\n", OUT_FILE, N*N*N);

    // ── Verification: spot checks ─────────────────────────────────────────────
    // Compare Bz written to the TNtuple against GetFieldXYZ called directly.
    // Differences should be at the float-vs-double rounding level (~1e-7 T).
    printf("=== Verification: spot checks ===\n");
    printf("  %-22s  %10s  %10s  %10s\n",
           "point (x,y,z) cm", "Bz direct", "Bz file", "delta");

    struct Spot { float x, y, z; };
    const Spot spots[] = {
        {  0.f,  0.f,   0.f },   // on-axis, z = 0
        {  0.f,  0.f, -24.f },   // on-axis, near measured Bz peak
        {  0.f,  0.f,   4.f },   // on-axis, near OPERA Bz peak
        { 30.f,  0.f,   0.f },   // mid-radius
        { 60.f,  0.f,   0.f },   // large radius
        {  0.f,  0.f, 100.f },   // z near edge of tracking volume
    };
    const int NSPOT = 6;

    // Compute expected Bz directly
    float bz_direct[NSPOT];
    for (int s = 0; s < NSPOT; ++s) {
        double Bx, By, Bz;
        meas.GetFieldXYZ(spots[s].x*CM2MM, spots[s].y*CM2MM, spots[s].z*CM2MM,
                         Bx, By, Bz);
        bz_direct[s] = (float)Bz;
    }

    // Re-open and find matching entries by direct index lookup.
    // Grid is x slowest/z fastest so entry = ix*N*N + iy*N + iz.
    TFile *fChk = TFile::Open(OUT_FILE, "READ");
    TNtuple *ntc = (TNtuple*)fChk->Get("fieldmap");
    float bzb;
    ntc->SetBranchAddress("bz", &bzb);

    for (int s = 0; s < NSPOT; ++s) {
        int ix = (int)std::lround((spots[s].x - GMIN) / DG);
        int iy = (int)std::lround((spots[s].y - GMIN) / DG);
        int iz = (int)std::lround((spots[s].z - GMIN) / DG);
        ntc->GetEntry(ix*N*N + iy*N + iz);
        printf("  (%+5.0f, %+5.0f, %+5.0f) cm    "
               "%+9.5f T  %+9.5f T  %+.2e T\n",
               spots[s].x, spots[s].y, spots[s].z,
               bz_direct[s], bzb, bz_direct[s] - bzb);
    }
    fChk->Close();

    printf("\nDone.\n");
}
