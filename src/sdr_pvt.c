//
//  Pocket SDR C Library - GNSS SDR PVT Functions
//
//  Author:
//  T.TAKASU
//
//  References:
//  [1] RINEX: The Receiver Independent Exchange Format version 3.05,
//      December 1, 2020
//
//  History:
//  2024-04-28  1.0  new
//  2024-12-30  1.1  add and update log contents
//                   add nav data consistency tests
//
#include "pocket_sdr.h"
#ifndef WIN32
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#endif

// constants and macros --------------------------------------------------------
#define SDR_EPOCH      1.0      // epoch time interval (s)
#define LAG_EPOCH      0.5      // max PVT epoch lag (s)
#define EL_MASK        15.0     // elavation mask (deg)
#define MAX_NOBS       256      // max number of obs data in a epoch
#define STD_ERR        0.015    // std-dev of carrier phase noise (m)
#define FILE_NAV       ".pocket_navdata.csv" // navigation data file
#define MAXDTGLO       (86400.0*7.0) // max age of GLONASS almanac (s)
#define MAXDTPOS       86400.0   // max age of last fix position (s)
#define SC_SPP_FAIL_RESET 10     // PPP KF reset after this many consecutive SPP failures in SC mode

#define ROUND(x)   (int)floor((x) + 0.5)
#define SQRT(x)    ((x) > 0.0 ? sqrt(x) : 0.0)
#define EQ(x, y)   (fabs((double)((x) - (y))) < 1e-12)

// global variable -------------------------------------------------------------
double sdr_epoch     = SDR_EPOCH;
double sdr_lag_epoch = LAG_EPOCH;
double sdr_el_mask   = EL_MASK;
double sdr_maxgdop   = 30.0;   /* reject threshold of gdop (default: 30.0) */
int    sdr_ionoopt   = IONOOPT_BRDC;
int    sdr_pmode     = PMODE_SINGLE;
double sdr_fixpos[3] = {0.0};  // known ECEF position for fixed-position mode (m)
int    sdr_ps_prn      = 0;     // pseudo-satellite PRN for AOWR time-transfer (0: disabled)
double sdr_ps_dist     = 0.6;   // true distance to pseudo-sat transmitter (m)
double sdr_ps_freq_err = 0.0;   // PS transmitter LO frequency error (Hz), e.g. -7.15e-3 for bladeRF
int    sdr_ps_gs_mode  = 0;     // GS mode: write clock_diff for SC (0: disabled)
char   sdr_ps_gs_file[256] = "./clock_diff.txt"; // output file for GS→SC IPC

// GS text ring-buffer — matches write_clock_difference() in gnss-sdr and
// getClockDiff() in jaxa-asyncOWR-prototype/src/bladeGPS/gpssim.c
// Format: 30 lines × 36 bytes, each line "tow(16),clock_diff_s(18)\n"
#define GS_LINE_SIZE  36
#define GS_NUM_LINES  30
#define GS_FILE_SIZE  (GS_LINE_SIZE * GS_NUM_LINES)   // 1080 bytes

// SC shared-data struct (binary) — matches hybrid_shared_data in gnss-sdr rtklib_pvt_gs
// Reserved for future SC-side implementation; not written in GS mode.
typedef struct {
    uint64_t seq;           // incremented each write; SC checks before/after for consistency
    double   tag_tow;       // estimated GPS TOW of PS transmission (s)
    double   clock_diff_s;  // PS clock offset from GPS = -dt_aowr_cp + rx_clock_offset (s)
    double   range_m;       // GS-to-PS distance (m), = sdr_ps_dist
} aowr_gs_shared_t;

static int   s_gs_fd       = -1;
static char *s_gs_map      = NULL;
static int   s_gs_cur_line = 0;
static double           s_rx_clock_s       = 0.0; // last receiver clock offset (s), from PVT epoch
static int              s_rx_sol_valid     = 0;   // 1 when PVT solution is valid
static double           s_gs_last_dt_aowr  = 0.0; // last valid dt_aowr_cp (s), frozen after PS loss
static int              s_gs_ps_ever_obs   = 0;   // 1 once PS has been observed (never reset)
int    sdr_ps_sc_mode  = 0;     // SC mode: read AOWR clock offset from file (0: disabled)
char   sdr_ps_sc_file[256] = "./dt_aowr_gnss.txt"; // input file: dt_aowr_gnss.txt from SC AOWR process

// SC binary seqlock file — matches hybrid_shared_data in gnss-sdr rtklib_pvt_gs.h
// Layout: 2 slots of aowr_gs_shared_t (64 bytes total); reader uses slot 0 only.
#define SC_FILE_SIZE       ((int)(sizeof(aowr_gs_shared_t) * 2))  // 64 bytes
static uint64_t s_sc_last_seq = 0; // last consumed seqlock sequence number
// Clock history for WLS drift estimation (mirrors linear_regression_by_wls in gnss-sdr).
#define SC_CLOCK_HIST_SIZE 3

static int    s_sc_fd          = -1;
static void  *s_sc_map         = NULL;
static double s_sc_last_dt_aowr = 0.0; // last valid dt_aowr_cp on SC side
static int    s_sc_ps_ever_obs  = 0;   // 1 once PS has been observed on SC side
static double s_sc_rx_clock     = 0.0; // latest AOWR-derived clock offset (s)
static double s_sc_last_tag_tow = 0.0; // tag_tow from the last AOWR file read
static double s_sc_clk_hist[SC_CLOCK_HIST_SIZE]; // rx clock offset at each reading
static double s_sc_tow_hist[SC_CLOCK_HIST_SIZE]; // SC GNSS TOW at each reading
static int    s_sc_hist_n      = 0;    // valid entries (0..SC_CLOCK_HIST_SIZE)
static double s_sc_dt_i        = 0.0;  // WLS intercept at tow_hist[0] (s)
static double s_sc_clock_drift = 0.0;  // WLS clock drift rate (s/s)
static rtk_t *s_sc_ref_rtk    = NULL; // reference PPP solver (uncorrected obs)
static sol_t  s_sc_ref_sol     = {0};  // last reference solution (SPP mode; mirror of s_sc_ref_rtk->sol in PPP mode)
static int    s_spp_fail_n     = 0;    // consecutive SPP_SEED failures in SC clock-fixed mode
static int    s_ref_spp_fail_n = 0;    // consecutive REF_SPP_SEED failures for s_sc_ref_rtk

/* SC antenna attitude — defined in lib/RTKLIB/src/rtkcmn.c (keeps librtk.a self-contained).
 * Set at startup from the options file; satazel() reads them in the antenna-frame path. */
extern int    sdr_sc_ant_fix;   // 0: standard ENU (terrestrial), 1: fixed SC antenna attitude
extern double sdr_sc_ant_az;    // boresight azimuth (deg)
extern double sdr_sc_ant_el;    // boresight elevation (deg; -90 = nadir toward Earth center)

int    sdr_dynamics  = 0;      // enable dynamics model for kinematic PPP
double sdr_prnaccelh = -1.0;   // horizontal acceleration process noise (m/s²), <0 = use RTKLIB default
double sdr_prnaccv   = -1.0;   // vertical acceleration process noise (m/s²), <0 = use RTKLIB default
int    sdr_dopvel    = 0;      // enable Doppler-based velocity update in kinematic PPP (0:off,1:on)
static const int systems[] = {
    SYS_GPS, SYS_GLO, SYS_GAL, SYS_QZS, SYS_CMP, SYS_IRN, SYS_SBS, 0
};

// satellite ID to system ------------------------------------------------------
static int sat2sys(const char *sat)
{
    static const char *str = "GREJCIS";
    const char *p = strchr(str, sat[0]);
    return p ? systems[(int)(p - str)] : 0;
}

// system to system index ------------------------------------------------------
static int sys2idx(int sys)
{
    for (int i = 0; systems[i]; i++) {
        if (systems[i] == sys) return i;
    }
    return -1;
}

// satellite to system index ---------------------------------------------------
static int sys_idx(int sat)
{
    return sys2idx(satsys(sat, NULL));
}

// signal ID to signal code -----------------------------------------------------
static uint8_t sig2code(const char *sig)
{
    static const char *sigs[] = {
        "L1CA" , "L1S"  , "L1CB" , "L1CP" , "L1CD" , "L2CM" , "L2CL" ,
        "L5I"  , "L5Q"  , "L5SI" , "L5SQ" , "L5SIV", "L5SQV", "L6D"  ,
        "L6E"  , "G1CA" , "G2CA" , "G1OCD", "G1OCP", "G2OCP", "G3OCD",
        "G3OCP", "E1B"  , "E1C"  , "E5AI" , "E5AQ" , "E5ABQ", "E5BI" ,
        "E5BQ" , "E6B"  , "E6C"  , "B1I"  , "B1CD" , "B1CP" , "B2I"  ,
        "B2AD" , "B2AP" , "B2BI" , "B3I"  , "I1SD" , "I1SP" , "I5S"  ,
        "ISS"  , NULL
    };
    static const uint8_t codes[] = {
        CODE_L1C, CODE_L1Z, CODE_L1E, CODE_L1L, CODE_L1S, CODE_L2S, CODE_L2L,
        CODE_L5I, CODE_L5Q, CODE_L5D, CODE_L5P, CODE_L5D, CODE_L5P, CODE_L6S,
        CODE_L6E, CODE_L1C, CODE_L2C, CODE_L4A, CODE_L4B, CODE_L6B, CODE_L3I,
        CODE_L3Q, CODE_L1B, CODE_L1C, CODE_L5I, CODE_L5Q, CODE_L8Q, CODE_L7I,
        CODE_L7Q, CODE_L6B, CODE_L6C, CODE_L2I, CODE_L1D, CODE_L1P, CODE_L7I,
        CODE_L5D, CODE_L5P, CODE_L7D, CODE_L6I, CODE_L1D, CODE_L1P, CODE_L5A,
        CODE_L9A
    };
    for (int i = 0; sigs[i]; i++) {
        if (!strcmp(sig, sigs[i])) return codes[i];
    }
    return 0;
}

//------------------------------------------------------------------------------
//  Output log $CH (receiver channel information).
//
//  format:
//      $CH,time,ch,rfch,sat,sig,prn,lock,cn0,coff,dop,adr,ssync,bsync,fsync,
//          rev,srev,err_phas,err_code,tow_v,tow,week,type,nnav,nerr,nlol,nfec
//          time  receiver time (s)
//          ch    receiver channel number
//          rfch  RF channel number
//          sat   satellite ID
//          sig   signal ID
//          prn   PRN number
//          lock  lock time (s)
//          cn0   C/N0 (dB-Hz)
//          coff  code offset (ms)
//          dop   Doppler frequency (Hz)
//          adr   accumlated Doppler range (cyc)
//          ssync secondary code sync flag (0:async, 1:sync)
//          bsync symbol/bit sync flag (0:async, 1:sync)
//          fsync frame sync flag (0:async, 1:sync)
//          rev   primary code polarity (0:normal, 1:reversed)
//          srev  secondary code polarity (0:normal, 1:reversed)
//          err_phas phase error (cyc)
//          err_code code error (10^-6 s)
//          tow_v  tow valid flag (0:invalid, 1:valid, 2:ambiguity unresolved)
//          tow   time of week (ms)
//          week  week number (week)
//          type  navigation subframe or message type
//          nnav  navigation subframe/message count
//          nerr  error subframe/message count
//          nlol  loss-of-lock count
//          nfec  number of error corrected (bits)
//
static void out_log_ch(sdr_ch_t *ch)
{
    sdr_log(3, "$CH,%.3f,%d,%d,%s,%s,%d,%.3f,%.1f,%.9f,%.3f,%.3f,%d,%d,%d,%d,"
        "%d,%.3f,%.3f,%d,%d,%d,%d,%d,%d,%d,%d", ch->time, ch->no, ch->rf_ch + 1,
        ch->sat, ch->sig, ch->prn, ch->lock * ch->T, ch->cn0, ch->coff * 1e3,
        ch->fd, ch->adr, ch->trk->sec_sync != 0, ch->nav->ssync != 0,
        ch->nav->fsync != 0, ch->nav->rev, ch->trk->sec_pol == -1,
        ch->trk->err_phas, ch->trk->err_code * 1e6, ch->tow_v, ch->tow,
        ch->week, ch->nav->type, ch->nav->count[0], ch->nav->count[1], ch->lost,
        ch->nav->nerr);
}

//------------------------------------------------------------------------------
//  Output log $OBS (observation data).
//
//  format:
//      $OBS,time,year,month,day,hour,min,sec,sat,code,cn0,pr,cp,dop,lli,fcn,ch
//          time  receiver time (s)
//          year,month,day  obs data day (GPST)
//          hour,min,sec  obs data time (GPST)
//          sat   satellite ID
//          code  RINEX obs code
//          cn0   C/N0 (dB-Hz)
//          pr    pseudorange (m)
//          cp    carrier phase (cyc)
//          dop   Doppler frequency (Hz)
//          lli   loss of lock indicator
//          fcn   frequency channel number for GLONASS FDMA signals
//          ch    RF ch number for L1 data
//
static void out_log_obs(double time, const obs_t *obs, const nav_t *nav)
{
    for (int i = 0; i < obs->n; i++) {
        const obsd_t *data = obs->data + i;
        double ep[6];
        char sat[16];
        time2epoch(data->time, ep);
        satno2id(data->sat, sat);
        if (sat[0] == '1') sat[0] = 'S';
        int fcn = 0, prn;
        if (satsys(data->sat, &prn) == SYS_GLO) {
            fcn = nav->geph[prn-1].frq;
        }
        for (int j = 0; j < NFREQ + NEXOBS; j++) {
            if (!data->code[j]) continue;
            sdr_log(3, "$OBS,%.3f,%.0f,%.0f,%.0f,%.0f,%.0f,%.3f,%s,%s,%.1f,"
                "%.3f,%.3f,%.3f,%d,%d,%d", time, ep[0], ep[1], ep[2], ep[3],
                ep[4], ep[5], sat, code2obs(data->code[j]),
                data->SNR[j] * SNR_UNIT, data->P[j], data->L[j], data->D[j],
                data->LLI[j], fcn, data->rcv);
        }
    }
}

//------------------------------------------------------------------------------
/* DOP for 3-unknown position-only solver (clock pre-corrected, not estimated).
   Returns dop[0]=0 (no GDOP), dop[1]=PDOP3, dop[2]=HDOP3, dop[3]=VDOP3.
   The H matrix has 3 columns (ENU: east, north, up) with no clock row, so
   there is no clock-position coupling. PDOP3 < PDOP (4-unknown) for the
   same geometry because the clock ambiguity does not inflate position error. */
static void dops_pos(int ns, const double *azel, double elmin, double *dop)
{
    double H[MAXOBS * 3], A[9] = {0};
    int i, r, c, n = 0;

    dop[0] = dop[1] = dop[2] = dop[3] = 0.0;
    for (i = 0; i < ns; i++) {
        if (azel[1 + i*2] < elmin) continue;
        double az = azel[i*2], el = azel[1 + i*2];
        H[n*3+0] = cos(el)*sin(az);   /* east */
        H[n*3+1] = cos(el)*cos(az);   /* north */
        H[n*3+2] = sin(el);           /* up */
        n++;
    }
    if (n < 3) return;
    /* A = H^T H (3×3, column-major: A[r + c*3]) */
    for (i = 0; i < n; i++)
        for (r = 0; r < 3; r++)
            for (c = 0; c < 3; c++)
                A[r + c*3] += H[i*3+r] * H[i*3+c];
    if (matinv(A, 3)) return;  /* singular → leave dop zero */
    /* dop[0]=0: no GDOP (clock not a state); dop[1]=PDOP3, dop[2]=HDOP3, dop[3]=VDOP3 */
    dop[1] = SQRT(A[0] + A[4] + A[8]); /* trace = east+north+up */
    dop[2] = SQRT(A[0] + A[4]);        /* horizontal */
    dop[3] = SQRT(A[8]);               /* vertical */
}

//  Output log $POS (position solution).
//
//  format:
//      $POS,time,year,month,day,hour,min,sec,lat,lon,hgt,Q,ns,stdn,stde,stdu,
//        dtr,x,y,z,vx,vy,vz,gdop,pdop,hdop,vdop
//          time  receiver time (s)
//          year,month,day  solution day (GPST)
//          hour,min,sec  solution time (GPST)
//          lat   solution latitude (deg, +:north, -:south)
//          lon   solution longitude (deg, +:east, -:west)
//          hgt   solution ellipsoidal height (m)
//          Q     quality flag (=5: single, =6: PPP)
//          ns    number of valid satellites
//          stdn  solution standard deviation north (m)
//          stde  solution standard deviation east (m)
//          stdu  solution standard deviation up (m)
//          dtr   receiver clock bias (s)
//          x,y,z ECEF position (m)
//          vx,vy,vz ECEF velocity (m/s, non-zero only in kinematic mode)
//          gdop  GDOP (0 when clock is pre-corrected — 3-unknown solve)
//          pdop  PDOP (3-unknown position-only when clock pre-corrected)
//          hdop  HDOP
//          vdop  VDOP
//          dtr_drift  clock drift rate (s/s): from Doppler in SPP; WLS in SC AOWR; 0 in PPP
//
static void out_log_pos(double time, const sol_t *sol, int nsat,
                        const ssat_t *ssat)
{
    double ep[6], pos[3], P[9], Q[9], dop[4] = {0};
    /* sol->dtr[0] unit depends on which solver produced the solution:
     *   pppos (stat=SOLQ_PPP)  → meters (= x[IC_GPS] from KF state)
     *   pntpos (stat=SOLQ_SINGLE) → seconds, even when pmode=PPP (SPP fallback)
     * PPP always uses uncorrected observations and estimates its own clock; dtr
     * is always in meters regardless of AOWR mode. */
    int dtr_in_meters = (sol->stat == SOLQ_PPP);
    double dtr_s      = dtr_in_meters ? sol->dtr[0] / CLIGHT : sol->dtr[0];
    time2epoch(timeadd(sol->time, dtr_s), ep);
    ecef2pos(sol->rr, pos);
    P[0] = sol->qr[0];
    P[4] = sol->qr[1];
    P[8] = sol->qr[2];
    P[1] = P[3] = sol->qr[3];
    P[5] = P[7] = sol->qr[4];
    P[2] = P[6] = sol->qr[5];
    covenu(pos, P, Q);
    if (ssat) {
        double azels[MAXSAT * 2];
        int ns = 0;
        for (int i = 0; i < MAXSAT; i++) {
            if (ssat[i].azel[1] > 0.0) {
                azels[ns * 2]     = ssat[i].azel[0];
                azels[ns * 2 + 1] = ssat[i].azel[1];
                ns++;
            }
        }
        /* When clock is pre-corrected (SC AOWR mode), the estimation has only 3
           position unknowns — use position-only DOP (no clock-position coupling).
           gdop=0 in the log signals "clock-fixed" mode to downstream tools.
           Otherwise compute standard 4-unknown GDOP/PDOP/HDOP/VDOP. */
        if (sdr_ps_sc_mode && s_sc_hist_n > 0)
            dops_pos(ns, azels, 0.0, dop);
        else
            dops(ns, azels, 0.0, dop);
    }
    sdr_log(3, "$POS,%.3f,%.0f,%.0f,%.0f,%.0f,%.0f,%.3f,%.9f,%.9f,%.3f,%d,%d,"
        "%.3f,%.3f,%.3f,%.9f,%.4f,%.4f,%.4f,%.6f,%.6f,%.6f,%.1f,%.1f,%.1f,%.1f,%.3e",
        time, ep[0], ep[1], ep[2], ep[3], ep[4], ep[5],
        pos[0] * R2D, pos[1] * R2D, pos[2], sol->stat, sol->ns, SQRT(Q[4]),
        SQRT(Q[0]), SQRT(Q[8]), dtr_s,
        sol->rr[0], sol->rr[1], sol->rr[2],
        sol->rr[3], sol->rr[4], sol->rr[5],
        dop[0], dop[1], dop[2], dop[3], sol->dtr[5]);
}

//------------------------------------------------------------------------------
//  Output log $ATT (attitude solution).
//
//  format:
//      $ATT,time,year,month,day,hour,min,sec,roll,pitch,yaw,nobs,bias1,bias2,
//          bias3,bias4,bias5,bias6,bias7,bias8
//          time  receiver time (s)
//          year,month,day  solution day (GPST)
//          hour,min,sec  solution time (GPST)
//          roll  roll angle (deg)
//          pitch pitch angle (deg)
//          yaw   yaw angle (deg)
//          nobs  number of obs data
//          bias1 hardware bias RFCH 1 (ns)
//          bias2 hardware bias RFCH 2 (ns)
//          ...
//          bias8 hardware bias RFCH 8 (ns)
//
#if 0
static void out_log_att(double time, const sdr_att_t *att)
{
    double ep[6];
    time2epoch(att->time, ep);
    sdr_log(3, "$ATT,%.3f,%.0f,%.0f,%.0f,%.0f,%.0f,%.3f,%.3f,%.3f,%.3f,%d,%.3f,"
        "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f", time, ep[0], ep[1], ep[2], ep[3],
        ep[4], ep[5], att->roll * R2D, att->pitch * R2D, att->yaw * R2D,
        att->pitch * R2D, att->nobs, att->bias[0] * CLIGHT / 1e9,
        att->bias[1] * CLIGHT / 1e9, att->bias[2] * CLIGHT / 1e9,
        att->bias[3] * CLIGHT / 1e9, att->bias[4] * CLIGHT / 1e9,
        att->bias[5] * CLIGHT / 1e9, att->bias[6] * CLIGHT / 1e9,
        att->bias[7] * CLIGHT / 1e9);
}
#endif

//------------------------------------------------------------------------------
//  Output log $SAT (satellite information).
//
//  format:
//      $SAT,time,year,month,day,hour,min,sec,sat,pvt,obs,cn0,az,el,res
//          time  receiver time (s)
//          year,month,day  solution day (GPST)
//          hour,min,sec  solution time (GPST)
//          sat   satellite ID
//          pvt   PVT status (0: not used, 1: used)
//          obs   L1 obs data status (0: not available, 1: available)
//          cn0   L1 signal C/N0 (dB-Hz)
//          az    azimuth angle (deg)
//          el    elavation angle (deg)
//          res   L1 pseudorange residual (m)
//
static void out_log_sat(double time, int sat, const sol_t *sol,
    const ssat_t *ssat)
{
    double ep[6];
    char str[16];
    time2epoch(timeadd(sol->time, sol->dtr[0]), ep);
    satno2id(sat, str);
    if (str[0] == '1') str[0] = 'S';
    sdr_log(3, "$SAT,%.3f,%.0f,%.0f,%.0f,%.0f,%.0f,%.3f,%s,%d,%d,%.1f,%.1f,"
        "%.1f,%.3f", time, ep[0], ep[1], ep[2], ep[3], ep[4], ep[5], str,
        ssat->vs, ssat->snr[0] > 0, ssat->snr[0] * SNR_UNIT,
        ssat->azel[0] * R2D, ssat->azel[1] * R2D, ssat->resp[0]);
}

//------------------------------------------------------------------------------
//  Output log $EPH (decoded ephemeris).
//
//  format:
//      $EPH,time,sat,sig,IODE,IODC,SVA,SVH,Toe,Toc,Ttr,A,e,i0,OMEGA0,omega,M0,
//        delta-n,OMEGAdot,Idot,Crc,Crs,Cuc,Cus,Cic,Cis,Toes,Fit,Af0,Af1,Af2,
//        TGD,code,flag (GPS,Galileo,QZSS,BeiDou,NavIC)
//      $EPH,time,sat,sig,tb,fcn,SVH,SVA,age,Toe,Tof,pos-x,pos-y,pos-z,vel-x,
//        vel-y,vel-z,acc-x,acc-y,acc-z,tau-n,gamma-n,delta-tau-n (GLONASS)
//      
//          time  receiver time (s)
//          sat   satellite ID
//          sig   signal ID
//          ...   ephemeris parameters
//
static void out_log_eph(double time, const char *sat, const char *sig,
    const void *p)
{
    char buff[2048];
    
    if (sat[0] != 'R') {
        const eph_t *eph = (const eph_t *)p;
        snprintf(buff, sizeof(buff), "%d,%d,%d,%d,%d,%d,%d,%.14E,%.14E,%.14E,"
            "%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,"
            "%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%d,%d", eph->iode,
            eph->iodc, eph->sva, eph->svh, (int)eph->toe.time, (int)eph->toc.time,
            (int)eph->ttr.time, eph->A, eph->e, eph->i0, eph->OMG0, eph->omg,
            eph->M0, eph->deln, eph->OMGd, eph->idot, eph->crc, eph->crs,
            eph->cuc, eph->cus, eph->cic, eph->cis, eph->toes, eph->fit,
            eph->f0, eph->f1, eph->f2, eph->tgd[0], eph->code, eph->flag);
    } else {
        const geph_t *geph = (const geph_t *)p;
        snprintf(buff, sizeof(buff), "%d,%d,%d,%d,%d,%d,%d,%.14E,%.14E,%.14E,"
            "%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E", geph->iode,
            geph->frq, geph->svh, geph->sva, geph->age, (int)geph->toe.time,
            (int)geph->tof.time, geph->pos[0], geph->pos[1], geph->pos[2],
            geph->vel[0], geph->vel[1], geph->vel[2], geph->acc[0], geph->acc[1],
            geph->acc[2], geph->taun, geph->gamn,geph->dtaun);
    }
    sdr_log(3, "$EPH,%.3f,%s,%s,%s", time, sat, sig, buff);
}

//------------------------------------------------------------------------------
//  Output log $ALM (decoded almanac).
//
//  format:
//      $ALM,time,sat,sig,svh,week,A,e,i0,OMG0,omg,M0,OMGd,toas,f0,f1
//      (GPS,Galileo,QZSS,BeiDou)
//          time   receiver time (s)
//          sat    satellite ID
//          sig    signal ID
//          svh    SV health (0=ok)
//          week   almanac week
//          A      semi-major axis (m)
//          e      eccentricity
//          i0     inclination (rad)
//          OMG0   RAAN (rad)
//          omg    arg of perigee (rad)
//          M0     mean motion (rad)
//          OMGd   RAAN rate (rad/s)
//          toas   toa (s)
//          f0,f1  clock bias/drift
//     $ALM,time,sat,sig,svh,taun,lambda,di,eps,omg,tlambda,dT,dTd,frq
//     (GLONASS)
//          time   receiver time (s)
//          sat    satellite ID
//          sig    signal ID
//          svh    SV health (0=ok)
//          taun   clock correction (s)
//          lambda ascending node longitude (rad)
//          di     inclination correction (rad)
//          eps    eccentricity
//          omg    arg of perigee (rad)
//          tlambda ascending node time (s)
//          dT,dTd period correction (s/orbit,s/orbit^2)
//          frq    frequency channel number
//
static void out_log_alm(double time, const char *sig, const alm_t *alm)
{
    char id[8];

    if (alm->sat <= 0) return;
    satno2id(alm->sat, id);
    if (id[0] == 'R') { // GLONASS (modified Keplerian almanac)
        sdr_log(3, "$ALM,%.3f,%s,%s,%d,%.10E,%.10E,%.10E,%.10E,%.10E,%.10E,"
            "%.10E,%.10E,%d", time, id, sig, alm->svh, alm->glo.taun,
            alm->glo.lambda, alm->glo.di, alm->glo.eps, alm->glo.omg,
            alm->glo.tlambda, alm->glo.dT, alm->glo.dTd, alm->glo.frq);
    } else {
        sdr_log(3, "$ALM,%.3f,%s,%s,%d,%d,%.10E,%.10E,%.10E,%.10E,%.10E,%.10E,"
            "%.10E,%.10E,%.10E,%.10E", time, id, sig, alm->svh, alm->week,
            alm->A, alm->e, alm->i0, alm->OMG0, alm->omg, alm->M0, alm->OMGd,
            alm->toas, alm->f0, alm->f1);
    }
}

// output all valid almanac of a satellite system ------------------------------
static void out_log_alm_sys(double time, const char *sig, const nav_t *nav,
    int sys)
{
    int sat, prn;

    for (sat = 1; sat <= MAXSAT; sat++) {
        if (satsys(sat, &prn) != sys) continue;
        if (nav->alm[sat-1].sat == sat) out_log_alm(time, sig, nav->alm + sat - 1);
    }
}

// output NMEA RMC, GGA, GSA and GSV -------------------------------------------
static void out_nmea(const sol_t *sol, const ssat_t *ssat, stream_t *str)
{
    uint8_t buff[4096];
    int n = 0;
    if (!str) return;
    n += outnmea_rmc(buff + n, sol);
    n += outnmea_gga(buff + n, sol);
    n += outnmea_gsa(buff + n, sol, ssat);
    n += outnmea_gsv(buff + n, sol, ssat);
    sdr_str_write(str, buff, n);
}

// count number of signals -----------------------------------------------------
//   rcv_no: RF CH filter (0: all)
static int num_sigs(int idx, const obs_t *obs, int rcv_no)
{
    int nsig = 0, mask[MAXCODE] = {0};
    
    for (int i = 0; i < obs->n; i++) {
        obsd_t *data = obs->data + i;
        if (sys_idx(data->sat) != idx) continue;
        if (rcv_no > 0 && data->rcv != rcv_no) continue;
        for (int j = 0; j < NFREQ + NEXOBS; j++) {
            if (!data->code[j] || mask[data->code[j]-1]) continue;
            mask[data->code[j]-1] = 1;
            nsig++;
        }
    }
    return nsig;
}

// output RTCM3 MSM messages for a single RF CH (rcv_no = 0: all CHs merged) ---
static void out_rtcm3_msm(rtcm_t *rtcm, const obs_t *obs, stream_t *str,
    int rcv_no)
{
    // RTCM3 MSM message types
    static const int msgs[] = {1077, 1087, 1097, 1117, 1127, 1137, 1107, 0};
    int nsig[7] = {0}, idx_tail = -1;
    
    for (int i = 0; msgs[i]; i++) {
        if ((nsig[i] = num_sigs(i, obs, rcv_no))) idx_tail = i;
    }
    if (idx_tail < 0) return;
    
    rtcm->staid = rcv_no; // 0:all,1-:rcv_no
    
    for (int i = 0; msgs[i]; i++) {
        if (!nsig[i]) continue;
        rtcm->obs.n = 0;
        for (int j = 0; j < obs->n; j++) {
            obsd_t *data = obs->data + j;
            if (sys_idx(data->sat) != i) continue;
            if (rcv_no > 0 && data->rcv != rcv_no) continue;
            
            // separate messages if nsat x nsig > 64
            if ((rtcm->obs.n + 1) * nsig[i] > 64) {
                if (gen_rtcm3(rtcm, msgs[i], 0, 1)) {
                    sdr_str_write(str, rtcm->buff, rtcm->nbyte);
                }
                rtcm->obs.n = 0;
            }
            rtcm->obs.data[rtcm->obs.n++] = *data;
        }
        if (rtcm->obs.n > 0 && gen_rtcm3(rtcm, msgs[i], 0, i < idx_tail)) {
            sdr_str_write(str, rtcm->buff, rtcm->nbyte);
        }
    }
}

// output RTCM3 observation data -----------------------------------------------
static void out_rtcm3_obs(rtcm_t *rtcm, const obs_t *obs, stream_t *str,
    const sdr_rcv_t *rcv)
{
    if (!str || obs->n <= 0) return;
    
    rtcm->time = obs->data[0].time;
    
    if (strstr(rcv->opt, "-ARRAY")) { // set staid by RF CH number
        int nch = rcv->nrfch + rcv->narch;
        for (int rcv_no = 1; rcv_no <= nch; rcv_no++) {
            out_rtcm3_msm(rtcm, obs, str, rcv_no);
        }
    } else {
        out_rtcm3_msm(rtcm, obs, str, 0);
    }
}

// output RTCM3 navigation data ------------------------------------------------
static void out_rtcm3_nav(rtcm_t *rtcm, int sat, int type, const nav_t *nav,
    stream_t *str)
{
    // RTCM3 navigation message types
    static const int msgs[] = {1019, 1020, 1046, 1044, 1042, 1041, 0, 0};
    int prn, sys = satsys(sat, &prn), idx = sys_idx(sat);
    
    if (!str || idx < 0 || !msgs[idx]) return;
    if (sys == SYS_GLO) {
        rtcm->nav.geph[prn-1] = nav->geph[prn-1];
    } else {
        rtcm->nav.eph[MAXSAT*type+sat-1] = nav->eph[MAXSAT*type+sat-1];
    }
    rtcm->ephsat = sat;
    int msg = (sys == SYS_GAL && type == 1) ? 1045 : msgs[idx];
    if (gen_rtcm3(rtcm, msg, 0, 0)) {
        sdr_str_write(str, rtcm->buff, rtcm->nbyte);
    }
}

// set observation data index --------------------------------------------------
static void set_obs_idx(sdr_rcv_t *rcv)
{
    int codes[7][NFREQ+NEXOBS] = {{0}};
    
    for (int i = 0; i < rcv->nch; i++) {
        sdr_ch_t *ch = rcv->th[i]->ch;
        int sys = sat2sys(ch->sat), code = sig2code(ch->sig);
        int j = sys2idx(sys), k = code2idx(sys, code);
        if (j < 0 || k < 0 || codes[j][k] == code) continue;
        if (codes[j][k] == 0) {
            codes[j][k] = code;
            continue;
        }
        for (k = NFREQ; k < NFREQ + NEXOBS; k++) {
            if (codes[j][k] == code) break;
            if (codes[j][k] == 0) {
                codes[j][k] = code;
                break;
            }
        }
    }
    for (int i = 0; i < rcv->nch; i++) {
        sdr_ch_t *ch = rcv->th[i]->ch;
        int sys = sat2sys(ch->sat), code = sig2code(ch->sig);
        int j = sys2idx(sys);
        if (j < 0) continue;
        for (int k = 0; k < NFREQ + NEXOBS; k++) {
            if (codes[j][k] != code) continue;
            ch->obs_idx = k;
            break;
        }
        if (sys == SYS_QZS && !strcmp(ch->sig, "L1CB")) {
            ch->obs_idx = 0; // L1C/A and L1C/B exclusive
        }
    }
}

// save nav data + almanac + last fix to file -----------------------------------
static void save_navdata(const char *file, const nav_t *nav,
    const sdr_pvt_t *pvt)
{
    FILE *fp;
    char id[8];
    int sat, prn;

    savenav(file, nav); // eph/geph/IONUTC (overwrites with "w")

    if (!(fp = fopen(file, "a"))) return; // append ALM and POS

    // almanac (all systems)
    for (sat = 1; sat <= MAXSAT; sat++) {
        if (nav->alm[sat-1].sat != sat) continue;
        satno2id(sat, id);
        if (id[0] == 'R') { // GLONASS
            fprintf(fp, "ALM,%s,%d,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,"
                "%.14E,%.14E,%d\n", id, nav->alm[sat-1].svh,
                nav->alm[sat-1].glo.taun, nav->alm[sat-1].glo.lambda,
                nav->alm[sat-1].glo.di,   nav->alm[sat-1].glo.eps,
                nav->alm[sat-1].glo.omg,  nav->alm[sat-1].glo.tlambda,
                nav->alm[sat-1].glo.dT,   nav->alm[sat-1].glo.dTd,
                nav->alm[sat-1].glo.frq);
        } else {
            if (satsys(sat, &prn) == SYS_SBS) continue;
            fprintf(fp, "ALM,%s,%d,%d,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,"
                "%.14E,%.14E,%.14E,%.14E\n", id, nav->alm[sat-1].svh,
                nav->alm[sat-1].week, nav->alm[sat-1].A,
                nav->alm[sat-1].e,    nav->alm[sat-1].i0,
                nav->alm[sat-1].OMG0, nav->alm[sat-1].omg,
                nav->alm[sat-1].M0,   nav->alm[sat-1].OMGd,
                nav->alm[sat-1].toas, nav->alm[sat-1].f0,
                nav->alm[sat-1].f1);
        }
    }
    // last fix position
    if (pvt->last_valid) {
        fprintf(fp, "POS,%ld,%.9f,%.4f,%.4f,%.4f,%.14E,%.14E\n",
            (long)pvt->last_time.time, pvt->last_time.sec,
            pvt->last_rr[0], pvt->last_rr[1], pvt->last_rr[2],
            pvt->last_dtr, pvt->last_dtrd);
    }
    fclose(fp);
}

// resolve almanac week from current CPU time ----------------------------------
static int resolve_alm_week(double toas)
{
    int week;
    double tow = time2gpst(utc2gpst(timeget()), &week);

    if      (toas < tow - 302400.0) week++;
    else if (toas > tow + 302400.0) week--;
    return week;
}

// test almanac week consistency with current CPU time -------------------------
static int valid_alm_week(int week)
{
    int w;
    (void)time2gpst(utc2gpst(timeget()), &w);
    return week > 0 && abs(week - w) <= 128;
}

// normalize GLONASS frequency channel number ---------------------------------
static int norm_glo_fcn(int fcn)
{
    return fcn > 15 ? fcn - 32 : fcn;
}

// load nav data + almanac + last fix from file ---------------------------------
static void load_navdata(const char *file, nav_t *nav, sdr_pvt_t *pvt)
{
    FILE *fp;
    char buff[4096], id[16];
    int sat, svh, week, frq;

    readnav(file, nav); // eph/geph/IONUTC

    if (!(fp = fopen(file, "r"))) return;

    while (fgets(buff, sizeof(buff), fp)) {
        if (!strncmp(buff, "ALM,", 4)) {
            if (sscanf(buff + 4, "%15[^,]", id) < 1) continue;
            if (!(sat = satid2no(id))) continue;
            if (sat < 1 || sat > MAXSAT) continue;
            nav->alm[sat-1].sat = sat;
            if (id[0] == 'R') { // GLONASS
                sscanf(buff + 4, "%*[^,],%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d",
                    &svh,
                    &nav->alm[sat-1].glo.taun,   &nav->alm[sat-1].glo.lambda,
                    &nav->alm[sat-1].glo.di,      &nav->alm[sat-1].glo.eps,
                    &nav->alm[sat-1].glo.omg,     &nav->alm[sat-1].glo.tlambda,
                    &nav->alm[sat-1].glo.dT,      &nav->alm[sat-1].glo.dTd,
                    &frq);
                nav->alm[sat-1].svh     = svh;
                nav->alm[sat-1].glo.frq = norm_glo_fcn(frq);
            } else {
                sscanf(buff + 4, "%*[^,],%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                    &svh, &week,
                    &nav->alm[sat-1].A,    &nav->alm[sat-1].e,
                    &nav->alm[sat-1].i0,   &nav->alm[sat-1].OMG0,
                    &nav->alm[sat-1].omg,  &nav->alm[sat-1].M0,
                    &nav->alm[sat-1].OMGd, &nav->alm[sat-1].toas,
                    &nav->alm[sat-1].f0,   &nav->alm[sat-1].f1);
                nav->alm[sat-1].svh  = svh;
                nav->alm[sat-1].week = valid_alm_week(week) ? week :
                    resolve_alm_week(nav->alm[sat-1].toas);
                nav->alm[sat-1].toa = gpst2time(nav->alm[sat-1].week,
                    nav->alm[sat-1].toas);
            }
        } else if (!strncmp(buff, "POS,", 4)) {
            long t;
            double sec, rr[3], dtr, dtrd;
            gtime_t time;
            if (sscanf(buff + 4, "%ld,%lf,%lf,%lf,%lf,%lf,%lf", &t, &sec,
                rr, rr+1, rr+2, &dtr, &dtrd) == 7) {
                time.time = t;
                time.sec  = sec;
                if (fabs(timediff(utc2gpst(timeget()), time)) < MAXDTPOS) {
                    pvt->last_time  = time;
                    pvt->last_rr[0] = rr[0];
                    pvt->last_rr[1] = rr[1];
                    pvt->last_rr[2] = rr[2];
                    pvt->last_dtr   = dtr;
                    pvt->last_dtrd  = dtrd;
                    pvt->last_valid = 1;
                } else {
                    pvt->last_valid = 0;
                }
            }
        }
    }
    fclose(fp);
}

// GS mode: open shared-data file and mmap it ---------------------------------
static void open_aowr_gs(void)
{
#ifndef WIN32
    s_gs_fd = open(sdr_ps_gs_file, O_RDWR | O_CREAT | O_TRUNC, 0644);
    if (s_gs_fd < 0) {
        sdr_log(1, "$LOG,0,,0,GS: cannot open %s: %s", sdr_ps_gs_file, strerror(errno));
        return;
    }
    // Extend file to GS_FILE_SIZE (matching gnss-sdr write_clock_difference init)
    if (lseek(s_gs_fd, GS_FILE_SIZE - 1, SEEK_SET) == -1 ||
        write(s_gs_fd, "", 1) == -1) {
        sdr_log(1, "$LOG,0,,0,GS: write error on %s", sdr_ps_gs_file);
        close(s_gs_fd); s_gs_fd = -1; return;
    }
    s_gs_map = (char *)mmap(NULL, GS_FILE_SIZE,
        PROT_READ | PROT_WRITE, MAP_SHARED, s_gs_fd, 0);
    if (s_gs_map == MAP_FAILED) {
        sdr_log(1, "$LOG,0,,0,GS: mmap failed on %s", sdr_ps_gs_file);
        s_gs_map = NULL; close(s_gs_fd); s_gs_fd = -1; return;
    }
    memset(s_gs_map, ' ', GS_FILE_SIZE);
    s_gs_cur_line = 0;
    sdr_log(3, "$LOG,0,,0,GS mode: sharing %s", sdr_ps_gs_file);
#endif
}

static void close_aowr_gs(void)
{
#ifndef WIN32
    if (s_gs_map) { munmap(s_gs_map, GS_FILE_SIZE); s_gs_map = NULL; }
    if (s_gs_fd >= 0) { close(s_gs_fd); s_gs_fd = -1; }
#endif
}

static void write_aowr_gs(double tag_tow, double clock_diff_s)
{
#ifndef WIN32
    if (!s_gs_map) return;
    // Text ring-buffer matching write_clock_difference() in gnss-sdr and
    // read by getClockDiff() in jaxa-asyncOWR-prototype:
    //   16-char tow + "," + 18-char clock_diff_s + "\n" = 36 bytes/line
    char tow_buf[32], dt_buf[32];
    snprintf(tow_buf, sizeof(tow_buf), "%16.15g", tag_tow);
    snprintf(dt_buf,  sizeof(dt_buf),  "%18.16g", clock_diff_s);
    tow_buf[16] = '\0';  // truncate to exactly 16 chars
    dt_buf[18]  = '\0';  // truncate to exactly 18 chars
    char line[GS_LINE_SIZE + 1];
    snprintf(line, sizeof(line), "%.16s,%.18s\n", tow_buf, dt_buf);
    memcpy(s_gs_map + (size_t)s_gs_cur_line * GS_LINE_SIZE, line, GS_LINE_SIZE);
    s_gs_cur_line = (s_gs_cur_line + 1) % GS_NUM_LINES;
#endif
}

// SC mode: open dt_aowr_gnss.txt for reading (O_RDONLY; file is written by SC AOWR process)
static void open_aowr_sc(void)
{
#ifndef WIN32
    s_sc_fd = open(sdr_ps_sc_file, O_RDONLY);
    if (s_sc_fd < 0) {
        sdr_log(1, "$LOG,0,,0,SC: cannot open %s: %s", sdr_ps_sc_file, strerror(errno));
        return;
    }
    s_sc_map = mmap(NULL, SC_FILE_SIZE, PROT_READ, MAP_SHARED, s_sc_fd, 0);
    if (s_sc_map == MAP_FAILED) {
        sdr_log(1, "$LOG,0,,0,SC: mmap failed on %s", sdr_ps_sc_file);
        s_sc_map = NULL; close(s_sc_fd); s_sc_fd = -1; return;
    }
    s_sc_hist_n = 0;
    /* Prime s_sc_last_seq with the file's current seq so stale data from a
     * previous run is skipped; only new GS writes after this point are read. */
    aowr_gs_shared_t *p0 = (aowr_gs_shared_t *)s_sc_map;
    s_sc_last_seq = p0->seq;
    sdr_log(3, "$LOG,0,,0,SC mode: reading AOWR clock from %s (hist=%d, skip_seq=%llu)",
        sdr_ps_sc_file, SC_CLOCK_HIST_SIZE, (unsigned long long)s_sc_last_seq);
#endif
}

static void close_aowr_sc(void)
{
#ifndef WIN32
    if (s_sc_map && s_sc_map != MAP_FAILED) {
        munmap(s_sc_map, SC_FILE_SIZE); s_sc_map = NULL;
    }
    if (s_sc_fd >= 0) { close(s_sc_fd); s_sc_fd = -1; }
#endif
}

// Seqlock read — mirrors read_hybrid_shared_data() in rtklib_pvt_gs.cc.
// Returns 1 with new data, 0 if write in progress or no new slot since last call.
static int read_aowr_sc(double *tag_tow, double *clock_diff_s, double *range_m)
{
#ifndef WIN32
    if (!s_sc_map) return 0;
    aowr_gs_shared_t *p = (aowr_gs_shared_t *)s_sc_map;
    uint64_t seq_before = p->seq;
    __sync_synchronize();
    aowr_gs_shared_t snap = *p;
    __sync_synchronize();
    uint64_t seq_after = p->seq;
    if (seq_before != seq_after) return 0;        // write in progress
    if (seq_after == 0 || seq_after == s_sc_last_seq) return 0;  // uninitialized or stale
    s_sc_last_seq = seq_after;
    *tag_tow      = snap.tag_tow;
    *clock_diff_s = snap.clock_diff_s;
    *range_m      = snap.range_m;
    return 1;
#else
    return 0;
#endif
}

//------------------------------------------------------------------------------
//  Generate a new SDR PVT.
//
//  args:
//      rcv      (I)  SDR receiver
//
//  returns:
//      SDR PVT (NULL: error)
//
sdr_pvt_t *sdr_pvt_new(sdr_rcv_t *rcv)
{
    sdr_pvt_t *pvt = (sdr_pvt_t *)sdr_malloc(sizeof(sdr_pvt_t));
    pvt->obs = (obs_t *)sdr_malloc(sizeof(obs_t));
    pvt->obs->data = (obsd_t *)sdr_malloc(sizeof(obsd_t) * MAX_NOBS);
    pvt->obs->nmax = MAX_NOBS;
    pvt->nav = (nav_t *)sdr_malloc(sizeof(nav_t));
    pvt->nav->eph = (eph_t *)sdr_malloc(sizeof(eph_t) * MAXSAT * 4);
    pvt->nav->n = pvt->nav->nmax = MAXSAT * 4;
    pvt->nav->geph = (geph_t *)sdr_malloc(sizeof(geph_t) * MAXPRNGLO);
    pvt->nav->ng = pvt->nav->ngmax = MAXPRNGLO;
    pvt->nav->seph = (seph_t *)sdr_malloc(sizeof(seph_t) * NSATSBS * 2);
    pvt->nav->ns = pvt->nav->nsmax = NSATSBS * 2;
    pvt->nav->alm = (alm_t *)sdr_malloc(sizeof(alm_t) * MAXSAT);
    pvt->nav->na = pvt->nav->namax = MAXSAT;
    pvt->sol = (sol_t *)sdr_malloc(sizeof(sol_t));
    pvt->ssat = (ssat_t *)sdr_malloc(sizeof(ssat_t) * MAXSAT);
    pvt->rtcm = (rtcm_t *)sdr_malloc(sizeof(rtcm_t));
    init_rtcm(pvt->rtcm);
    if (sdr_pmode >= PMODE_PPP_KINEMA) {
        prcopt_t opt = prcopt_default;
        opt.mode    = sdr_pmode;
        opt.navsys  = SYS_GPS | SYS_GLO | SYS_GAL | SYS_QZS | SYS_CMP | SYS_IRN;
        opt.ionoopt = (sdr_ionoopt==IONOOPT_GRAPHIC) ? IONOOPT_GRAPHIC : IONOOPT_IFLC;
        opt.tropopt = TROPOPT_ESTG;
        opt.maxinno = 100.0; /* relax outlier gate: default 30m too tight during PPP convergence */
        if (sdr_ionoopt==IONOOPT_GRAPHIC) {
            opt.modear = ARMODE_OFF; /* float bias only */
        }
        opt.elmin   = sdr_el_mask * D2R;
        opt.maxgdop = sdr_maxgdop;
        opt.dynamics = sdr_dynamics;
        if (sdr_prnaccelh >= 0.0) opt.prn[3] = sdr_prnaccelh;
        if (sdr_prnaccv   >= 0.0) opt.prn[4] = sdr_prnaccv;
        opt.dopvel = sdr_dopvel;
        if (sdr_pmode == PMODE_PPP_FIXED) {
            opt.ru[0] = sdr_fixpos[0];
            opt.ru[1] = sdr_fixpos[1];
            opt.ru[2] = sdr_fixpos[2];
        }
        pvt->rtk = (rtk_t *)sdr_malloc(sizeof(rtk_t));
        rtkinit(pvt->rtk, &opt);
        /* Reference PPP solver for SC AOWR evaluation (same opts, independent KF state) */
        if (sdr_ps_sc_mode) {
            s_sc_ref_rtk = (rtk_t *)sdr_malloc(sizeof(rtk_t));
            rtkinit(s_sc_ref_rtk, &opt);
        }
    }
    set_obs_idx(rcv);
    pvt->rcv = rcv;
    sdr_mutex_init(&pvt->mtx);
    load_navdata(FILE_NAV, pvt->nav, pvt); // load nav + almanac + last fix
    if (sdr_ps_gs_mode) open_aowr_gs();
    if (sdr_ps_sc_mode) open_aowr_sc();
    return pvt;
}

//------------------------------------------------------------------------------
//  Free a SDR PVT.
//
// Load precise navigation file (SP3 orbit or RINEX CLK) into pvt->nav. -------
// Called from pocket_trk.c after receiver open, before data arrives.
// Multiple calls are allowed: load SP3 first, then CLK for higher clock rate.
//
// Supported extensions (case-insensitive):
//   .sp3 .eph → precise orbit (+ coarse 15-min clock embedded in SP3)
//   .clk .clk_05s .clk_15s .clk_30s .clk_30 → RINEX clock (high-rate)
//
// When an SP3 file is loaded, rtk->opt.sateph is set to EPHOPT_PREC so
// pppos() interpolates the SP3 orbit instead of using broadcast ephemeris.
//
void sdr_pvt_loadnav(sdr_pvt_t *pvt, const char *file)
{
    if (!pvt || !file || !*file) return;

    // lowercase copy of the filename for extension matching
    char lc[2048];
    int k = 0;
    for (const char *p = file; *p && k < (int)sizeof(lc) - 1; p++)
        lc[k++] = tolower((unsigned char)*p);
    lc[k] = '\0';

    int is_sp3 = (strstr(lc, ".sp3") || strstr(lc, ".eph")) ? 1 : 0;
    int is_clk = strstr(lc, ".clk") ? 1 : 0;

    if (is_sp3) {
        int ne_before = pvt->nav->ne;
        readsp3(file, pvt->nav, 0);
        int ne_loaded = pvt->nav->ne - ne_before;
        if (ne_loaded > 0) {
            if (pvt->rtk) pvt->rtk->opt.sateph = EPHOPT_PREC;
            if (s_sc_ref_rtk) s_sc_ref_rtk->opt.sateph = EPHOPT_PREC;
            sdr_log(3, "$LOG,0.000,SP3 LOADED: %s (%d orbit epochs, sateph=PREC)",
                file, ne_loaded);
        } else {
            fprintf(stderr, "SP3 load failed or empty: %s\n", file);
        }
    } else if (is_clk) {
        int nc_before = pvt->nav->nc;
        readrnxc(file, pvt->nav);
        int nc_loaded = pvt->nav->nc - nc_before;
        if (nc_loaded > 0) {
            sdr_log(3, "$LOG,0.000,CLK LOADED: %s (%d clock epochs)", file, nc_loaded);
        } else {
            fprintf(stderr, "CLK load failed or empty: %s\n", file);
        }
    } else {
        fprintf(stderr, "pocket_trk -nav: unrecognised file type: %s\n", file);
    }
}

//  args:
//      pvt      (I)  SDR PVT generated by sdr_pvt_new()
//
//  returns:
//      none
//
void sdr_pvt_free(sdr_pvt_t *pvt)
{
    if (!pvt) return;
    close_aowr_gs();
    close_aowr_sc();
    save_navdata(FILE_NAV, pvt->nav, pvt); // save nav + almanac + last fix
    sdr_free(pvt->obs->data);
    sdr_free(pvt->obs);
    sdr_free(pvt->nav->eph);
    sdr_free(pvt->nav->geph);
    sdr_free(pvt->nav->seph);
    sdr_free(pvt->nav);
    sdr_free(pvt->sol);
    sdr_free(pvt->ssat);
    free_rtcm(pvt->rtcm);
    sdr_free(pvt->rtcm);
    if (pvt->rtk) {
        rtkfree(pvt->rtk);
        sdr_free(pvt->rtk);
    }
    if (s_sc_ref_rtk) {
        rtkfree(s_sc_ref_rtk);
        sdr_free(s_sc_ref_rtk);
        s_sc_ref_rtk = NULL;
    }
    sdr_free(pvt);
}

// initialize epoch time and cycle ---------------------------------------------
static void init_epoch(sdr_pvt_t *pvt, int64_t ix, sdr_ch_t *ch)
{
    if (!ch->week) return;
    // Never seed the receiver's absolute epoch reference from the PS/AOWR
    // ranging channel -- its week (when it decodes one at all) is not
    // trustworthy (see gen_prng()'s use_rx_week). In practice a real
    // satellite always decodes its week first (PS takes ~200s+ vs. <30s for
    // GPS), but this closes the gap regardless of timing.
    if (sdr_ps_prn > 0) {
        char ps_sat_id[16];
        sdr_sat_id("L1CA", sdr_ps_prn, ps_sat_id);
        if (!strcmp(ch->sat, ps_sat_id) && !strcmp(ch->sig, "L1CA")) return;
    }
    double tow = floor(ch->tow * 1e-3 / sdr_epoch) * sdr_epoch + sdr_epoch;
    pvt->time = gpst2time(ch->week, tow);
    pvt->ix = ix + ROUND((tow - ch->tow * 1e-3 - 0.07) / SDR_CYC);
    pvt->ix = (pvt->ix / 20) * 20; // round by 20 ms
}

// generate pseudorange --------------------------------------------------------
// use_rx_week: for the PS/AOWR ranging channel, ch->week decoded from that
// signal's own nav data is not trustworthy (see update_aowr()'s week-fold
// comment) -- always reconstruct tau from the receiver's own independently-
// known GPS week instead, regardless of whether this channel ever decodes
// (or mis-decodes) a week number of its own. Real satellite channels pass 0
// and keep the normal ch->week-trusting behavior.
static double gen_prng(gtime_t time, const sdr_ch_t *ch, int use_rx_week)
{
    int week;
    double tau = 0.0, tow = time2gpst(time, &week);

    if (use_rx_week) {
        tau = tow - ch->tow * 1e-3 + ch->coff;
        if (tau < -302400.0) tau += 604800.0;
        if (tau >  302400.0) tau -= 604800.0;
    } else if (ch->week > 0) {
        tau = (week - ch->week) * 86400.0 * 7 + tow - ch->tow * 1e-3 + ch->coff;
    } else if (ch->tow_v == 1) { // tow valid but GPS week not yet decoded from nav
        // use current GPS week from receiver time; handle end-of-week wrap
        tau = tow - ch->tow * 1e-3 + ch->coff;
        if (tau < -302400.0) tau += 604800.0;
        if (tau >  302400.0) tau -= 604800.0;
    } else if (ch->tow_v == 2) { // resolve 100 ms ambiguity (0.05 <= tau < 0.15)
        tau = tow - ch->tow * 1e-3 + ch->coff + ch->nav->coff;
        tau -= floor(tau / 0.1) * 0.1;
        if (tau < 0.05) tau += 0.1;
    }
    // for debug
    trace(3, "%s %-5s %3d %4d %10.3f %10.3f %12.9f %12.9f\n", ch->sat, ch->sig,
        ch->prn, ch->week, tow, ch->tow * 1e-3, ch->coff, tau);

    return CLIGHT * (tau + 0.5 * ch->T * ch->fd / ch->fc);
}

// generate carrier-phase ------------------------------------------------------
static double gen_cphas(const sdr_ch_t *ch, double P)
{
    double L = -ch->adr;
    
    L += (ch->nav->rev ? 0.5 : 0.0) + (ch->trk->sec_pol == 1 ? 0.5 : 0.0);
    
    // phase alignment ([1] Table A23)
    if (!strcmp(ch->sig, "L1CD") || !strcmp(ch->sig, "L1CP")) {
        L += 0.25; // + 1/4 cyc
    } else if (!strcmp(ch->sig, "L5Q") || !strcmp(ch->sig, "L5SQ") ||
        !strcmp(ch->sig, "L5SQV")) {
        L -= 0.25; // - 1/4 cyc
    } else if (!strcmp(ch->sig, "G3OCP") || !strcmp(ch->sig, "E5AQ") ||
        !strcmp(ch->sig, "E5ABQ") || !strcmp(ch->sig, "E5BQ") ||
        !strcmp(ch->sig, "B1CP") ||
        !strcmp(ch->sig, "B2AP")) {
        L += 0.25; // + 1/4 cyc
    } else if (!strcmp(ch->sig, "L2CM")) {
        L += (ch->sat[0] == 'J') ? 0.0 : -0.25; // 0 cyc (QZSS), -1/4 cyc (GPS)
    } else if ((!strcmp(ch->sig, "B1I") || !strcmp(ch->sig, "B2I")) &&
        (ch->prn <= 5 || ch->prn >= 59)) {
        L += 0.5;
    }
    return L;
}

// update observation data -----------------------------------------------------
static void update_obs(gtime_t time, obs_t *obs, sdr_ch_t *ch)
{
    uint8_t code = sig2code(ch->sig);
    double P = gen_prng(time, ch, 0);
    int i, idx = ch->obs_idx, sat;
    
    if (strstr(ch->sat, "R-") || strstr(ch->sat, "R+")) return;
    if (P <= 0.0 || idx < 0 || !(sat = satid2no(ch->sat))) return;
    
    for (i = 0; i < obs->n; i++) {
        if (sat == obs->data[i].sat &&
            (idx != 0 || obs->data[i].rcv == ch->rf_ch + 1)) break;
    }
    if (i >= obs->n) {
        if (i >= obs->nmax) return;
        memset(obs->data + i, 0, sizeof(obsd_t));
        obs->data[i].time = time;
        obs->data[i].sat = sat;
        obs->data[i].rcv = idx == 0 ? ch->rf_ch + 1 : 0;
        obs->n++;
    }
    obs->data[i].code[idx] = code;
    obs->data[i].P[idx] = P;
    obs->data[i].L[idx] = gen_cphas(ch, P);
    obs->data[i].D[idx] = (float)ch->fd;
    obs->data[i].SNR[idx] = (uint16_t)(ch->cn0 / SNR_UNIT + 0.5);

    // Correct carrier phase and Doppler for non-zero IF offset.
    // The PLL drives ch->fd = fd_satellite + fi (total baseband, not just Doppler)
    // because phi = fi*tau + adr uses only the last-interval fi, forcing the PLL
    // to fold cumulative fi into fd. Left uncorrected: GF drifts -0.91 m/s
    // (false cycle slips every epoch), MW drifts +3.2 m/s, and Lc innovations
    // grow 0.35 m/s — all preventing carrier-phase PPP.
    // Fix: remove accumulated fi component from carrier phase and Doppler.
    // Applies to any channel with non-zero IF (e.g. fi_L1=+1.831 Hz, fi_L2=-2.289 Hz
    // for PocketSDR FE with external 10 MHz reference).
    if (ch->fi != 0.0) {
        obs->data[i].L[idx] += ch->fi * (ch->lock * ch->T);  // remove fi*t cycles
        obs->data[i].D[idx] -= (float)ch->fi;                  // remove fi from Doppler
    }
    if (ch->lock * ch->T <= 2.0 || fabs(ch->trk->err_phas) > 0.25) {
        obs->data[i].LLI[idx] |= 1; // PLL unlock
    }
    if (ch->nav->fsync <= 0 && ch->trk->sec_sync <= 0) {
        obs->data[i].LLI[idx] |= 2; // half-cyc-amb unknown
    }
}

static void update_aowr(double time, gtime_t gtime, double P, double L);

//------------------------------------------------------------------------------
//  Update observation data.
//
//  args:
//      pvt      (IO) SDR PVT
//      ix       (I)  received IF data cycle (cyc)
//      ch       (IO) SDR receiver channel
//
//  returns:
//      none
//
void sdr_pvt_udobs(sdr_pvt_t *pvt, int64_t ix, sdr_ch_t *ch)
{
    sdr_mutex_lock(&pvt->mtx);

    if (pvt->ix <= 0) { // initialize epoch time and cycle
        init_epoch(pvt, ix, ch);
    }
    if (pvt->ix > 0 && ix == pvt->ix) { // update observation data
        if (ch->state == SDR_STATE_LOCK && ch->tow >= 0 && ch->tow_v > 0 &&
            (ch->nav->fsync > 0 || ch->trk->sec_sync > 0)) {
            update_obs(pvt->time, pvt->obs, ch);
        }
        pvt->nch++;
    }
    // Output $CH log at sdr_epoch intervals.
    int64_t log_step = (int64_t)(sdr_epoch / SDR_CYC + 0.5);
    int log_now = (pvt->ix > 0) ?
        (ix == pvt->ix) : (log_step > 0 && ix % log_step == 0);
    
    if (log_now && ch->state == SDR_STATE_LOCK && ch->lock > 0) {
        out_log_ch(ch);
    }
    // 20 ms high-rate AOWR tick (50 Hz), independent of 1 Hz PVT epoch
    if (sdr_ps_prn > 0 && pvt->ix > 0 && ix % 20 == 0) {
        char ps_sat_id[16];
        sdr_sat_id("L1CA", sdr_ps_prn, ps_sat_id);
        if (!strcmp(ch->sat, ps_sat_id) && !strcmp(ch->sig, "L1CA") &&
            ch->state == SDR_STATE_LOCK && ch->tow >= 0 && ch->tow_v > 0 &&
            (ch->nav->fsync > 0 || ch->trk->sec_sync > 0)) {
            gtime_t gt = timeadd(pvt->time, (ix - pvt->ix) * SDR_CYC);
            double P = gen_prng(gt, ch, 1); // always use receiver's own week for PS
            double L = gen_cphas(ch, P);
            if (ch->fi != 0.0)         L += ch->fi         * (ch->lock * ch->T);
            if (sdr_ps_freq_err != 0.0) L += sdr_ps_freq_err * (ch->lock * ch->T);
            update_aowr(ix * SDR_CYC, gt, P, L);
        }
    }
    sdr_mutex_unlock(&pvt->mtx);
}

// test nav data consistency for GLONASS ---------------------------------------
static int test_nav_glo(const sdr_ch_t *ch)
{
    int t[3];
    for (int i = 0; i < 3; i++) {
        t[i] = ch->nav->lock_sf[i+1] - ch->nav->lock_sf[i];
    }
    return t[0] == 2000 && t[1] == 2000 && t[2] == 2000;
}

// test match of ephemeris parameters ------------------------------------------
static int match_eph(const eph_t *e1, const eph_t *e2)
{
    return EQ(e1->iode, e2->iode) && EQ(e1->iodc, e2->iodc) &&
        EQ(e1->A, e2->A) && EQ(e1->e, e2->e) && EQ(e1->i0, e2->i0) &&
        EQ(e1->OMG0, e2->OMG0) && EQ(e1->omg, e2->omg) && EQ(e1->M0, e2->M0) &&
        EQ(e1->deln, e2->deln) && EQ(e1->OMGd, e2->OMGd) &&
        EQ(e1->idot, e2->idot) && EQ(e1->crc, e2->crc) &&
        EQ(e1->crs, e2->crs) && EQ(e1->cuc, e2->cuc) && EQ(e1->cus, e2->cus) &&
        EQ(e1->cic, e2->cic) && EQ(e1->cis, e2->cis) && EQ(e1->f0, e2->f0) &&
        EQ(e1->f1, e2->f1) && EQ(e1->f2, e2->f2) &&
        EQ(e1->tgd[0], e2->tgd[0]) && EQ(e1->toes, e2->toes);
}

// test nav data consistency for BeiDou D1/D2 ----------------------------------
static int test_match_eph(eph_t *eph1, const eph_t *eph2)
{
    if (match_eph(eph1 + MAXSAT, eph2)) { // match previous ephemeris
        *(eph1 + MAXSAT) = *eph2;
        *eph1 = *eph2;
        return 1;   
    } else { // not match previous ephemeris
        *(eph1 + MAXSAT) = *eph2;
        return 0;   
    }
}

//------------------------------------------------------------------------------
//  Update navigation data.
//
//  args:
//      pvt      (IO) SDR PVT
//      ch       (IO) SDR receiver channel
//
//  returns:
//      none
//
void sdr_pvt_udnav(sdr_pvt_t *pvt, sdr_ch_t *ch)
{
    uint8_t *data = ch->nav->data;
    int prn, sat = satid2no(ch->sat), sys = satsys(sat, &prn);
    
    if (sys == SYS_NONE) return;
    
    sdr_mutex_lock(&pvt->mtx);
    
    if (!strcmp(ch->sig, "L1CA") && sys == SYS_SBS) { // SBAS
        if (ch->nav->type == 9) { // geo navigation message
            int week, tow = (int)time2gpst(pvt->time, &week);
            sbsmsg_t msg = {week, tow, (uint8_t)ch->prn, 1};
            memcpy(msg.msg, data, 29);
            if (sbsupdatecorr(&msg, pvt->nav) == 9) {
                pvt->count[2]++;
            }
        }
    } else if (!strcmp(ch->sig, "L1CA") || !strcmp(ch->sig, "L1CB")) { // GPS/QZS LNAV
        if (ch->nav->type == 3 &&
            decode_frame(data, pvt->nav->eph + sat - 1, NULL, NULL, NULL)) {
            pvt->nav->eph[sat-1].sat = sat;
            out_log_eph(ch->time, ch->sat, ch->sig, pvt->nav->eph + sat - 1);
            out_rtcm3_nav(pvt->rtcm, sat, 0, pvt->nav, pvt->rcv->strs[1]);
            pvt->count[2]++;
        }
        if (sys == SYS_GPS && ch->nav->type == 4) {
            decode_frame(data, NULL, NULL, pvt->nav->ion_gps, NULL);
        }
        if ((ch->nav->type == 4 || ch->nav->type == 5) && // almanac (SF4/SF5)
            decode_frame(data, NULL, pvt->nav->alm, NULL, NULL)) {
            out_log_alm_sys(ch->time, ch->sig, pvt->nav,
                sys == SYS_QZS ? SYS_QZS : SYS_GPS);
        }
    } else if (!strcmp(ch->sig, "G1CA") || !strcmp(ch->sig, "G2CA")) { // GLO NAV
        pvt->nav->geph[prn-1].tof = pvt->time;
        if (ch->nav->type == 4 && test_nav_glo(ch) &&
            decode_glostr(data, pvt->nav->geph + prn - 1, NULL, NULL)) {
            pvt->nav->geph[prn-1].sat = sat;
            pvt->nav->geph[prn-1].frq = ch->prn; // FCN
            out_log_eph(ch->time, ch->sat, ch->sig, pvt->nav->geph + prn - 1);
            out_rtcm3_nav(pvt->rtcm, sat, 0, pvt->nav, pvt->rcv->strs[1]);
            pvt->count[2]++;
        }
        if (ch->nav->type == 15 && // almanac (strings 6-15 complete)
            decode_glostr(data, NULL, pvt->nav->alm, NULL)) {
            out_log_alm_sys(ch->time, ch->sig, pvt->nav, SYS_GLO);
        }
    } else if (!strcmp(ch->sig, "E1B") || !strcmp(ch->sig, "E5BI")) { // GAL I/NAV
        if (ch->nav->type == 4 &&
            decode_gal_inav(data, pvt->nav->eph + sat - 1, NULL, NULL, NULL)) {
            pvt->nav->eph[sat-1].sat = sat;
            out_log_eph(ch->time, ch->sat, ch->sig, pvt->nav->eph + sat - 1);
            out_rtcm3_nav(pvt->rtcm, sat, 0, pvt->nav, pvt->rcv->strs[1]);
            pvt->count[2]++;
        }
        if (ch->nav->type == 10 && // almanac (word types 7-10 complete)
            decode_gal_inav(data, NULL, pvt->nav->alm, NULL, NULL)) {
            out_log_alm_sys(ch->time, ch->sig, pvt->nav, SYS_GAL);
        }
    } else if (!strcmp(ch->sig, "E5AI")) { // GAL F/NAV
        if (ch->nav->type == 4 &&
            decode_gal_fnav(data, pvt->nav->eph + MAXSAT + sat - 1, NULL,
                NULL, NULL)) {
            pvt->nav->eph[MAXSAT+sat-1].sat = sat;
            out_log_eph(ch->time, ch->sat, ch->sig, pvt->nav->eph + MAXSAT +
                sat - 1);
            out_rtcm3_nav(pvt->rtcm, sat, 1, pvt->nav, pvt->rcv->strs[1]);
            pvt->count[2]++;
        }
        if (ch->nav->type == 6 && // almanac (page types 5-6 complete)
            decode_gal_fnav(data, NULL, pvt->nav->alm, NULL, NULL)) {
            out_log_alm_sys(ch->time, ch->sig, pvt->nav, SYS_GAL);
        }
    } else if (!strcmp(ch->sig, "B1I") || !strcmp(ch->sig, "B2I") ||
             !strcmp(ch->sig, "B3I")) {
        eph_t eph = {0};
        if (ch->prn >= 6 && ch->prn <= 58) { // BDS D1 NAV
            if (ch->nav->type == 3 && decode_bds_d1(data, &eph, NULL, NULL, NULL)) {
                if (test_match_eph(pvt->nav->eph + sat - 1, &eph)) {
                    pvt->nav->eph[sat-1].sat = sat;
                    out_log_eph(ch->time, ch->sat, ch->sig, pvt->nav->eph + sat - 1);
                    out_rtcm3_nav(pvt->rtcm, sat, 0, pvt->nav, pvt->rcv->strs[1]);
                    pvt->count[2]++;
                } else {
                    out_log_eph(ch->time, ch->sat, ch->sig, &eph);
                    sdr_log(3, "$LOG,%.3f,%s,%s,EPHEMERIS UNMATCH", ch->time,
                        ch->sat, ch->sig);
                }
            }
            if ((ch->nav->type == 4 || ch->nav->type == 5) && // almanac (SF4/SF5)
                decode_bds_d1(data, NULL, pvt->nav->alm, NULL, NULL)) {
                out_log_alm_sys(ch->time, ch->sig, pvt->nav, SYS_CMP);
            }
        } else { // BDS D2 NAV
            if (ch->nav->type == 10 && decode_bds_d2(data, &eph, NULL, NULL)) {
                if (test_match_eph(pvt->nav->eph + sat - 1, &eph)) {
                    pvt->nav->eph[sat-1].sat = sat;
                    out_log_eph(ch->time, ch->sat, ch->sig, pvt->nav->eph + sat - 1);
                    out_rtcm3_nav(pvt->rtcm, sat, 0, pvt->nav, pvt->rcv->strs[1]);
                    pvt->count[2]++;
                } else {
                    out_log_eph(ch->time, ch->sat, ch->sig, &eph);
                    sdr_log(3, "$LOG,%.3f,%s,%s,EPHEMERIS UNMATCH", ch->time,
                        ch->sat, ch->sig);
                }
            }
            if (ch->nav->type >= 100 && // almanac (SF5 page)
                decode_bds_d2(data, NULL, pvt->nav->alm, NULL)) {
                out_log_alm_sys(ch->time, ch->sig, pvt->nav, SYS_CMP);
            }
        }
    } else if (!strcmp(ch->sig, "I5S") || !strcmp(ch->sig, "ISS")) { // NavIC NAV
        if (ch->nav->type == 2 &&
            decode_irn_nav(data, pvt->nav->eph + sat - 1, NULL, NULL)) {
            pvt->nav->eph[sat-1].sat = sat;
            out_log_eph(ch->time, ch->sat, ch->sig, pvt->nav->eph + sat - 1);
            out_rtcm3_nav(pvt->rtcm, sat, 0, pvt->nav, pvt->rcv->strs[1]);
            pvt->count[2]++;
        }
    }
    sdr_mutex_unlock(&pvt->mtx);
}

// correct solution time -------------------------------------------------------
static void corr_sol_time(sol_t *sol)
{
    if (fabs(sol->dtr[0]) >= 1e-9) return;
    
    // use GLOT, GALT, BDT or IRT as solution time in case of GPS absence
    for (int i = 1; i < 5; i++) {
        if (fabs(sol->dtr[i]) < 1e-9) continue;
        sol->dtr[0] = sol->dtr[i];
        sol->time = timeadd(sol->time, -sol->dtr[0]);
        return;
    }
}

// update satellite az/el angles -----------------------------------------------
static void update_azel(const nav_t *nav, const sol_t *sol, ssat_t *ssat)
{
    for (int i = 0; i < MAXSAT; i++) {
        double rs[6], dts[2], var, pos[3], e[3];
        int svh;
        
        if (satpos(sol->time, sol->time, i + 1, EPHOPT_BRDC, nav, rs, dts,
            &var, &svh) && geodist(rs, sol->rr, e) > 0.0) {
            ecef2pos(sol->rr, pos);
            satazel(pos, e, ssat[i].azel);
        }
    }
}

// save last fix and output solution logs --------------------------------------
static void output_sol(sdr_pvt_t *pvt, double time)
{
    corr_sol_time(pvt->sol);
    pvt->last_time  = pvt->sol->time;
    pvt->last_rr[0] = pvt->sol->rr[0];
    pvt->last_rr[1] = pvt->sol->rr[1];
    pvt->last_rr[2] = pvt->sol->rr[2];
    pvt->last_dtr   = pvt->sol->dtr[0];
    pvt->last_dtrd  = pvt->sol->dtr[5];
    pvt->last_valid = 1;
    out_log_pos(time, pvt->sol, pvt->obs->n, pvt->ssat);
    out_nmea(pvt->sol, pvt->ssat, pvt->rcv->strs[0]);
    pvt->count[0]++;
    for (int i = 0; i < MAXSAT; i++) {
        if (pvt->ssat[i].snr[0] == 0) continue;
        out_log_sat(time, i + 1, pvt->sol, pvt->ssat + i);
    }
}

static void res_obs_amb(obs_t *obs, int sys, uint8_t code, double sec);

//------------------------------------------------------------------------------
//  Update AOWR inter-system time bias from pseudo-satellite observation.
//  Implements the PR+CP time-transfer algorithm from rtklib_pvt_gs.cc:3069-3163.
//
//  format:
//      $AOWR,time,year,month,day,hour,min,sec,prn,P,L,dt_raw,dt_pr,dt_cp,
//          count,outlier
//          time    receiver time (s)
//          year,month,day,hour,min,sec  GPST epoch
//          prn     pseudo-satellite PRN
//          P       pseudorange (m)
//          L       carrier phase (cycles)
//          dt_raw  P/CLIGHT - raw propagation delay (s)
//          dt_pr   smoothed PR-based GNSSR-AOWR time offset (s)
//          dt_cp   carrier-phase enhanced GNSSR-AOWR time offset (s)
//          count   valid (non-outlier) sample count
//          outlier 1=outlier epoch, 0=valid epoch
//
static void update_aowr(double time, gtime_t gtime, double P, double L)
{
    static const double DT_DEV_THRESH   = 3.0 / CLIGHT; // 10 ns gate
    static const double T_WARMUP  = 20.0; // skip DLL settling transient (s)
    static const int    DEV_COUNT_THRESH = 100;
    static int     initialized   = 0;
    static int64_t dt_int_s      = 0;
    static double  dt_frac_sum   = 0.0;
    static double  dt0_frac_sum  = 0.0;
    static int     count         = 0;
    static double  dt_aowr       = 0.0;
    static double  dt_aowr_cp    = 0.0;
    static double  cp_thresh     = 3.0 / CLIGHT;
    static double  diff_total    = 0.0;
    static int     dev_count     = 0;
    static double  dt_new_frac_sum  = 0.0;
    static double  dt0_new_frac_sum = 0.0;
    static int     dt_new_count     = 0;
    static double  diff_new_total   = 0.0;
    static double  initial_rx    = 0.0;
    double dt_current  = P / CLIGHT;
    double Ci          = L / FREQ1;   // carrier phase in light-seconds

    if (!initialized) {
        dt_int_s    = (int64_t)round(dt_current);
        initial_rx  = time;
        initialized = 1;
    }

    /* Strip GPS week-number contribution from the PS pseudorange.
     * The PS navigation message may carry a wrong week number (e.g. stale
     * firmware value), causing P/c to jump by an integer multiple of 604800 s.
     * Folding dt_current to within ±302400 s of the initial reference makes
     * the computation depend only on TOW, not the week number. */
    {
        double delta = dt_current - (double)dt_int_s;
        delta -= 604800.0 * round(delta / 604800.0);
        dt_current = (double)dt_int_s + delta;
    }

    double dt0 = count > 0 ? (double)dt_int_s + dt0_frac_sum / count : 0.0;
    double dt0_current = dt_current - sdr_ps_dist / CLIGHT - Ci;

    /* The DLL's narrow (0.25 Hz) noise bandwidth gives it a ~10-12 s settling
     * time (2nd-order loop, zeta=0.707) before the code-phase estimate (and
     * thus P) converges after channel lock. Exclude this transient from the
     * dt_pr/dt_cp/dt0/cp_thresh statistics so it cannot bias the average. */
    int warmup = (time - initial_rx) < T_WARMUP;

    int outlier = !warmup && (dt_aowr != 0.0) &&
        (fabs(dt_current - dt_aowr)     > DT_DEV_THRESH ||
         fabs(dt0_current - dt0)        > DT_DEV_THRESH ||
         fabs(dt0 + Ci - dt_aowr_cp)    > cp_thresh     ||
         (time - initial_rx) > T_WARMUP + 40.0);

    if (warmup) {
        // skip stats accumulation during the DLL settling transient
    } else if (outlier) {
        // A fresh streak of consecutive outliers starts a fresh candidate cluster;
        // stale accumulators from a previous (interrupted or failed) streak must
        // not leak in, or dt_new_frac_sum/dt0_new_frac_sum end up summed over more
        // epochs than dev_count divides by, making the "cluster" average nonsense
        // and preventing dt_new_count from ever reaching DEV_COUNT_THRESH.
        if (dev_count == 0) {
            dt_new_frac_sum  = 0.0;
            dt0_new_frac_sum = 0.0;
            diff_new_total   = 0.0;
            dt_new_count     = 0;
        }
        dev_count++;
        dt_new_frac_sum  += dt_current  - (double)dt_int_s;
        dt0_new_frac_sum += dt0_current - (double)dt_int_s;
        double dt_new   = (double)dt_int_s + dt_new_frac_sum / dev_count;
        double diff_new = fabs(dt_current - dt_new);
        diff_new_total += diff_new;
        if (dt_new != 0.0 && diff_new < DT_DEV_THRESH)
            dt_new_count++;
        else
            dt_new_count = 0;
    } else {
        dev_count      = 0;
        dt_frac_sum   += dt_current - (double)dt_int_s;
        count++;
        dt_aowr        = (double)dt_int_s + dt_frac_sum / count;

        dt0_frac_sum  += dt0_current - (double)dt_int_s;
        dt0            = (double)dt_int_s + dt0_frac_sum / count;

        if (dt_aowr_cp != 0.0) {
            diff_total += fabs(dt0 + Ci - dt_aowr_cp);
            cp_thresh   = 3.0 * diff_total / count;
        }
        dt_aowr_cp = dt0 + Ci;
    }

    if (dev_count >= DEV_COUNT_THRESH) {
        if (dt_new_count >= DEV_COUNT_THRESH) {
            dt_frac_sum      = dt_new_frac_sum;
            dt0_frac_sum     = dt0_new_frac_sum; /* keep dt0 consistent with the new count */
            count            = dt_new_count;
            dt_aowr          = (double)dt_int_s + dt_new_frac_sum / dt_new_count;
            dt_new_count     = 0;
            /* dt_aowr_cp/cp_thresh/diff_total are carrier-smoothed quantities that
             * were never tracked for the candidate cluster (only the raw-pseudorange
             * dt_new_frac_sum was). Clear them so the next accepted (non-outlier)
             * epoch re-derives dt_aowr_cp = dt0+Ci from scratch against the new
             * dt0_frac_sum/count, instead of comparing against the stale value from
             * the cluster we just abandoned. */
            dt_aowr_cp       = 0.0;
            diff_total       = 0.0;
            cp_thresh        = 3.0 / CLIGHT;
        }
        dev_count        = 0;
        dt_new_frac_sum  = 0.0;
        dt0_new_frac_sum = 0.0;
        diff_new_total   = 0.0;
    }

    // dt_aowr_cp is momentarily 0.0 for the one epoch a reanchor commits on (it is
    // re-derived from dt0+Ci on the next accepted epoch) -- skip publishing that
    // placeholder so consumers never see a bogus zero clock offset.
    if (sdr_ps_gs_mode && count > 0 && dt_aowr_cp != 0.0) {
        s_gs_last_dt_aowr = dt_aowr_cp;
        s_gs_ps_ever_obs  = 1;
    }
    if (sdr_ps_sc_mode && count > 0 && dt_aowr_cp != 0.0) {
        s_sc_last_dt_aowr = dt_aowr_cp;
        s_sc_ps_ever_obs  = 1;
    }

    double ep[6];
    time2epoch(gtime, ep);
    sdr_log(3, "$AOWR,%.3f,%.0f,%.0f,%.0f,%.0f,%.0f,%.3f,%d,%.3f,%.3f,"
        "%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d",
        time, ep[0], ep[1], ep[2], ep[3], ep[4], ep[5],
        sdr_ps_prn, P, L, dt_current, dt_aowr, dt_aowr_cp, dt0, cp_thresh,
        count, outlier);
}

// update PVT solution ---------------------------------------------------------
static void update_sol(sdr_pvt_t *pvt)
{
    double time = pvt->ix * SDR_CYC;
    obsd_t obs[MAXSAT];
    obsd_t obs_ref[MAXSAT]; /* uncorrected obs for SC reference solver */
    int mask[MAXSAT] = {0}, nobs = 0;
    int nobs_ref = 0, run_ref_sol = 0, sc_clk_ready = 0;
    double clock_est = 0.0; /* AOWR-derived clock applied to obs (0 when not in SC mode) */
    char msg[128] = "";

    // deduplicate: one entry per satellite (L1+L2 merged in same obsd_t)
    for (int i = 0; i < pvt->obs->n && nobs < MAXSAT; i++) {
        int sat = pvt->obs->data[i].sat;
        if (pvt->obs->data[i].P[0] == 0.0 || mask[sat-1]) continue;
        obs[nobs++] = pvt->obs->data[i];
        mask[sat-1] = 1;
    }
    // merge L2+ from split entries (race: L2CM thread runs before L1CA thread)
    for (int i = 0; i < pvt->obs->n; i++) {
        const obsd_t *d = pvt->obs->data + i;
        if (d->P[0] != 0.0) continue; // L1 present → already merged above
        for (int j = 0; j < nobs; j++) {
            if (obs[j].sat != d->sat) continue;
            for (int k = 1; k < NFREQ + NEXOBS; k++) {
                if (d->code[k] && !obs[j].code[k]) {
                    obs[j].code[k] = d->code[k];
                    obs[j].P[k]    = d->P[k];
                    obs[j].L[k]    = d->L[k];
                    obs[j].D[k]    = d->D[k];
                    obs[j].SNR[k]  = d->SNR[k];
                    obs[j].LLI[k]  = d->LLI[k];
                }
            }
            break;
        }
    }
    // Resolve msec pseudorange ambiguities after merging L1+L2 into one entry.
    // Must run here so L1CA reference (P[0]) is available in the same obsd_t when
    // res_obs_amb searches for it. Running before merge zeroes L2CM pseudorange
    // on epochs where the L2CM thread updates pvt->obs before the L1CA thread.
    {
        obs_t merged_obs;
        merged_obs.data = obs;
        merged_obs.n = nobs;
        merged_obs.nmax = MAXSAT;
        res_obs_amb(&merged_obs, SYS_GPS | SYS_QZS, CODE_L2S, 1e-3);
        res_obs_amb(&merged_obs, SYS_GPS | SYS_QZS, CODE_L5Q, 20e-3);
        res_obs_amb(&merged_obs, SYS_QZS, CODE_L5P, 20e-3);
        res_obs_amb(&merged_obs, SYS_GLO, CODE_L3Q, 10e-3);
        res_obs_amb(&merged_obs, SYS_SBS, CODE_L5Q, 2e-3);
        out_log_obs(time, &merged_obs, pvt->nav);
        out_rtcm3_obs(pvt->rtcm, &merged_obs, pvt->rcv->strs[1], pvt->rcv);
        // update_aowr is now called at 20 ms rate from sdr_pvt_udobs()
    }

    // SC mode: always save uncorrected obs for the reference solver.
    // Reference solver runs on raw obs regardless of AOWR clock availability.
    if (sdr_ps_sc_mode && nobs > 0) {
        memcpy(obs_ref, obs, sizeof(obsd_t) * nobs);
        nobs_ref    = nobs;
        run_ref_sol = 1;
    }

    // SC mode: extrapolate AOWR clock to current epoch using WLS drift estimation,
    // then pre-correct obs so pntpos/pppos solve position only (3 unknowns).
    // Mirrors linear_regression_by_wls + rx_clock_offset_est in rtklib_pvt_gs.cc.
    if (sdr_ps_sc_mode && s_sc_ps_ever_obs) {
        // (a) Check for new AOWR data and update clock history
        double tag_tow, clock_diff_s, range_m;
        if (read_aowr_sc(&tag_tow, &clock_diff_s, &range_m)) {
            double new_clock = s_sc_last_dt_aowr + clock_diff_s;
            /* SC GNSS TOW: tag_tow (PS tx time) + SC pseudorange offset.
             * May be negative due to GPS week rollover — wrap to [0, 604800). */
            double new_tow = tag_tow + s_sc_last_dt_aowr;
            while (new_tow <      0.0) new_tow += 604800.0;
            while (new_tow >= 604800.0) new_tow -= 604800.0;
            s_sc_last_tag_tow = tag_tow;

            sdr_log(3, "$LOG,%.3f,SC_AOWR_RAW tag_tow=%.3f clock_diff=%.9f"
                " last_dt=%.9f new_clk=%.9f",
                time, tag_tow, clock_diff_s, s_sc_last_dt_aowr, new_clock);

            /* Append to history; shift oldest out when full */
            if (s_sc_hist_n < SC_CLOCK_HIST_SIZE) {
                s_sc_tow_hist[s_sc_hist_n] = new_tow;
                s_sc_clk_hist[s_sc_hist_n] = new_clock;
                s_sc_hist_n++;
            } else {
                memmove(s_sc_tow_hist, s_sc_tow_hist+1,
                        (SC_CLOCK_HIST_SIZE-1)*sizeof(double));
                memmove(s_sc_clk_hist, s_sc_clk_hist+1,
                        (SC_CLOCK_HIST_SIZE-1)*sizeof(double));
                s_sc_tow_hist[SC_CLOCK_HIST_SIZE-1] = new_tow;
                s_sc_clk_hist[SC_CLOCK_HIST_SIZE-1] = new_clock;
            }
            s_sc_rx_clock = new_clock;

            /* (b) WLS linear regression: y = dt_i + drift*(t - tow_hist[0]) */
            if (s_sc_hist_n > 1) {
                double S0=0,S1=0,S2=0,T0=0,T1=0,t0=s_sc_tow_hist[0];
                for (int k=0; k<s_sc_hist_n; k++) {
                    double t = s_sc_tow_hist[k] - t0;
                    double y = s_sc_clk_hist[k];
                    S0+=1; S1+=t; S2+=t*t; T0+=y; T1+=t*y;
                }
                double det = S0*S2 - S1*S1;
                if (fabs(det) > 1e-12) {
                    s_sc_dt_i        = (S2*T0 - S1*T1) / det;
                    s_sc_clock_drift = (-S1*T0 + S0*T1) / det;
                }
            } else {
                s_sc_dt_i        = new_clock;
                s_sc_clock_drift = 0.0;
            }
        }

        /* (c) Extrapolate AOWR clock to current epoch using obs[0].time as GPS TOW. */
        if (s_sc_rx_clock != 0.0 && s_sc_hist_n > 0 && nobs > 0 &&
                s_sc_last_tag_tow > 0.0) {
            int wn;
            double current_tow = time2gpst(obs[0].time, &wn);
            double dt = current_tow - s_sc_tow_hist[0];
            if (dt < 0.0) dt += 604800.0;
            clock_est    = s_sc_dt_i + s_sc_clock_drift * dt;
            sc_clk_ready = 1;

            sdr_log(3, "$LOG,%.3f,SC_AOWR clk=%.9f drift=%.3e dt=%.3f hist=%d",
                time, clock_est, s_sc_clock_drift, dt, s_sc_hist_n);
        }
    }

    if (sdr_pmode >= PMODE_PPP_KINEMA) {
        // Set time interval for KF process noise propagation (rtkpos() normally does this)
        gtime_t obs_time = nobs > 0 ? obs[0].time : pvt->time;
        pvt->rtk->tt = pvt->rtk->sol.time.time > 0 ?
            timediff(obs_time, pvt->rtk->sol.time) : 0.0;

        sol_t spp_sol = {0};
        if (sdr_ps_sc_mode && !sc_clk_ready) {
            /* SC mode but no AOWR clock yet: suppress main solver, keep KF frozen */
            pvt->rtk->sol.stat = 0;
            pvt->rtk->sol.time = obs_time; /* advance time for next epoch's tt */
            sdr_log(3, "$LOG,%.3f,SPP_SEED SUPPRESSED (SC no clock)", time);
        } else {
        // Seed pppos EKF with SPP position (pppos needs rtk->sol.rr != {0,0,0})
        // This also sets rtk->sol.time = obs_time for next epoch's tt computation.
        prcopt_t spopt = prcopt_default;
        spopt.navsys |= SYS_GLO | SYS_GAL | SYS_QZS | SYS_CMP | SYS_IRN;
        spopt.elmin   = sdr_el_mask * D2R;
        spopt.maxgdop = sdr_maxgdop;
        // dtr[0..NSYS-1] holds meters when the previous epoch was PPP (SOLQ_PPP).
        // pntpos expects seconds; without this conversion all residuals are ~300 km off
        // and pntpos rejects every satellite as an outlier ("lack of valid sats").
        if (pvt->rtk->sol.stat == SOLQ_PPP) {
            for (int i = 0; i < NSYS; i++) pvt->rtk->sol.dtr[i] /= CLIGHT;
        }
        if (sc_clk_ready) {
            spopt.clock_bias_fixed = 1; /* P/L pre-corrected → 3-unknown solve */
            pvt->rtk->opt.clock_bias_fixed = 1; /* For PPP */
            pvt->rtk->sol.dtr[0] = clock_est; /* fixed clock bias (seconds) */
        }
        // Pass pvt->rtk->ssat so pntpos sets ssat[sat].vs=1; pppos needs vs=1 to accept obs
        // In fixed-position mode, seed pntpos from the known position for a better clock estimate.
        if (sdr_pmode == PMODE_PPP_FIXED && norm(sdr_fixpos, 3) > 1.0) {
            pvt->rtk->sol.rr[0] = sdr_fixpos[0];
            pvt->rtk->sol.rr[1] = sdr_fixpos[1];
            pvt->rtk->sol.rr[2] = sdr_fixpos[2];
        }
        // One-time bootstrap: seed position from REF only when SC position is
        // still zero (never solved yet). After first convergence the SC solver
        // is independent; it must NOT be overridden by REF every epoch.
        if (sc_clk_ready && s_sc_ref_rtk && norm(s_sc_ref_rtk->sol.rr, 3) > 1.0
            && norm(pvt->rtk->sol.rr, 3) < 1.0) {
            sdr_log(3, "$LOG,%.3f,SPP_SEED first seed from REF %.1f %.1f %.1f",
                time, s_sc_ref_rtk->sol.rr[0],
                s_sc_ref_rtk->sol.rr[1], s_sc_ref_rtk->sol.rr[2]);
            pvt->rtk->sol.rr[0] = s_sc_ref_rtk->sol.rr[0];
            pvt->rtk->sol.rr[1] = s_sc_ref_rtk->sol.rr[1];
            pvt->rtk->sol.rr[2] = s_sc_ref_rtk->sol.rr[2];
        }
        // Save a reliable position seed. pntpos may converge to a wrong position
        // when an outlier satellite dominates the WLS; use the KF position (which
        // is maintained by pppos independently) so the RAIM-FDE fallback below
        // always starts from a sane initial point.
        double spp_rr0[3];
        for (int j = 0; j < 3; j++)
            spp_rr0[j] = norm(pvt->rtk->x, 3) > 1.0 ? pvt->rtk->x[j]
                                                      : pvt->rtk->sol.rr[j];

        pntpos(obs, nobs, pvt->nav, &spopt, &pvt->rtk->sol, NULL, pvt->rtk->ssat, msg);
        if (pvt->rtk->sol.stat) spp_sol = pvt->rtk->sol;
        sdr_log(3, "$LOG,%.3f,SPP_SEED stat=%d dtr=%.9f msg=%s",
            time, pvt->rtk->sol.stat, pvt->rtk->sol.dtr[0],
            pvt->rtk->sol.stat ? "ok" : msg);

        // RAIM-FDE fallback: when all-satellite pntpos fails (e.g., a newly acquired
        // satellite with a code-period ambiguity causes WLS divergence), try excluding
        // each satellite in turn.  udclk_ppp() reinitialises the PPP KF clock from
        // sol.dtr[0] every epoch, so a wrong SPP clock contaminates pppos as well.
        // Excluding the outlier satellite restores a valid dtr and unblocks pppos.
        if (!pvt->rtk->sol.stat && nobs > 5) {
            for (int exc = 0; exc < nobs && !pvt->rtk->sol.stat; exc++) {
                obsd_t oe[MAXOBS]; int ne = 0;
                for (int j = 0; j < nobs; j++) if (j != exc) oe[ne++] = obs[j];
                sol_t st = pvt->rtk->sol;
                for (int j = 0; j < 3; j++) st.rr[j] = spp_rr0[j];
                pntpos(oe, ne, pvt->nav, &spopt, &st, NULL, pvt->rtk->ssat, msg);
                if (st.stat) {
                    pvt->rtk->sol = st;
                    sdr_log(3, "$LOG,%.3f,SPP_RAIM excl=%d ok dtr=%.9f",
                        time, exc, st.dtr[0]);
                }
            }
            if (pvt->rtk->sol.stat) spp_sol = pvt->rtk->sol;
        }
        /* PPP clock protection: when RAIM-FDE fails (two simultaneous outlier sats)
         * but the PPP KF clock is already converged, restore sol.dtr[0] from the KF
         * state before pppos() runs.  udclk_ppp() reinitialises x[IC] from
         * sol.dtr[0] every epoch; a wrong SPP dtr shifts x[IC] by tens of km,
         * causing all pre-fit residuals to exceed maxinno=30m and a permanent
         * PPPOS NO SOLUTION.  Only apply when not in clock_bias_fixed mode (which
         * sets dtr[0] from the external AOWR clock, not SPP). */
        if (!pvt->rtk->sol.stat && !sc_clk_ready) {
            int np_clk = pvt->rtk->opt.dynamics ? 9 : 3; /* IC(0) = np + 0 */
            double clk_var = pvt->rtk->P[np_clk + np_clk * pvt->rtk->nx];
            if (clk_var > 0.0 && clk_var < 100.0) { /* PPP converged: std < 10m */
                pvt->rtk->sol.dtr[0] = pvt->rtk->x[np_clk] / CLIGHT;
                sdr_log(3, "$LOG,%.3f,PPP_CLK_PROTECT dtr=%.9f (SPP/RAIM failed, clk_std=%.3f m)",
                    time, pvt->rtk->sol.dtr[0], sqrt(clk_var));
            }
        }

        /* In clock_bias_fixed mode, a failed pntpos corrupts sol.rr via diverged
         * Newton iterations.  Restore it to the pre-pntpos seed so the next epoch
         * starts from the same clean position instead of drifting further away. */
        if (!pvt->rtk->sol.stat && sc_clk_ready) {
            pvt->rtk->sol.rr[0] = spp_rr0[0];
            pvt->rtk->sol.rr[1] = spp_rr0[1];
            pvt->rtk->sol.rr[2] = spp_rr0[2];
        }

        /* PPP KF reset after consecutive SPP failures in SC clock-fixed mode.
         * When SPP fails (AOWR clock has drifted), pppos() keeps propagating the
         * position state unconstrained. After SC_SPP_FAIL_RESET epochs the KF
         * position and phase biases are too inconsistent to re-converge; zeroing
         * the full state lets pppos re-initialise cleanly on the next good epoch. */
        if (sc_clk_ready) {
            if (pvt->rtk->sol.stat) {
                s_spp_fail_n = 0;
            } else if (++s_spp_fail_n >= SC_SPP_FAIL_RESET) {
                int nx = pvt->rtk->nx;
                memset(pvt->rtk->x, 0, nx * sizeof(double));
                memset(pvt->rtk->P, 0, nx * nx * sizeof(double));
                s_spp_fail_n = 0;
                sdr_log(3, "$LOG,%.3f,PPP_KF_RESET (SPP failed %d epochs in SC mode)",
                    time, SC_SPP_FAIL_RESET);
            }
        }

        // Pre-seed clock states so P[IC] is non-zero before pppos.
        // udclk_ppp() unconditionally resets x[IC]=CLIGHT*sol.dtr[0] (white-noise
        // clock model), so x[IC] is always written inside pppos regardless.  This
        // block only matters for P[IC]: after PPP_KF_RESET, P is zeroed, and
        // udclk_ppp resets P[IC]=VAR_CLK.  The block is therefore redundant but
        // kept for clarity — it mirrors what udclk_ppp does, before pppos runs.
        {
            int np = pvt->rtk->opt.dynamics ? 9 : 3;
            for (int i = 0; i < NSYS; i++) {
                int ic = np + i;
                if (pvt->rtk->x[ic] == 0.0 && pvt->rtk->sol.dtr[i] != 0.0) {
                    double dtr = i == 0 ? pvt->rtk->sol.dtr[0] :
                                          pvt->rtk->sol.dtr[0] + pvt->rtk->sol.dtr[i];
                    pvt->rtk->x[ic] = CLIGHT * dtr;
                    pvt->rtk->P[ic + ic * pvt->rtk->nx] = 60.0 * 60.0; /* VAR_CLK */
                }
            }
        }

        // Carrier phase and Doppler are corrected for non-zero IF in update_obs()
        // (both L1CA and L2CM independently, relative to zero-IF reference).
        // GF/MW/Lc drift are all zero after correction. gf_drift stays 0.

        // PPP: dual-frequency Kalman filter via pppos(), but only when there are
        // enough dual-frequency observations. With ndual=0 there are no IFLC pairs
        // and pppos returns immediately, but the EKF position state was already
        // propagated forward with process noise → diverges. Skip entirely.
        {
            int ndual_pre = 0;
            for (int i = 0; i < nobs; i++) {
                if (obs[i].code[1] && obs[i].P[1] != 0.0 && obs[i].L[1] != 0.0) ndual_pre++;
            }
            if (ndual_pre >= 4) pppos(pvt->rtk, obs, nobs, pvt->nav);
            else sdr_log(3, "$LOG,%.3f,PPPOS SKIPPED (ndual=%d)", time, ndual_pre);
        }
        } /* end else (sc_clk_ready) */

        // Reference PPP: normal PPP on uncorrected obs for SC AOWR evaluation.
        // Independent KF state in s_sc_ref_rtk; result logged as $REFPOS.
        if (run_ref_sol && s_sc_ref_rtk) {
            gtime_t rtime = nobs_ref > 0 ? obs_ref[0].time : pvt->time;
            s_sc_ref_rtk->tt = s_sc_ref_rtk->sol.time.time > 0 ?
                timediff(rtime, s_sc_ref_rtk->sol.time) : 0.0;

            // Seed reference RTK position from normal SPP on uncorrected obs
            prcopt_t rspopt = prcopt_default;
            rspopt.navsys |= SYS_GLO | SYS_GAL | SYS_QZS | SYS_CMP | SYS_IRN;
            rspopt.elmin   = sdr_el_mask * D2R;
            rspopt.maxgdop = sdr_maxgdop;
            pntpos(obs_ref, nobs_ref, pvt->nav, &rspopt, &s_sc_ref_rtk->sol,
                   NULL, s_sc_ref_rtk->ssat, msg);
            sdr_log(3, "$LOG,%.3f,REF_SPP_SEED stat=%d dtr=%.9f",
                time, s_sc_ref_rtk->sol.stat, s_sc_ref_rtk->sol.dtr[0]);

            /* PPP KF reset after consecutive SPP failures, mirroring the main
             * solver's SC_SPP_FAIL_RESET logic above. Without this, a single bad
             * SPP fix (e.g. a marginal 4-satellite/high-PDOP solution) permanently
             * seeds udpos_ppp()'s one-time init with garbage, every subsequent
             * epoch's residuals get rejected as outliers, and s_sc_ref_rtk never
             * recovers for the rest of the session. */
            if (s_sc_ref_rtk->sol.stat) {
                s_ref_spp_fail_n = 0;
            } else if (++s_ref_spp_fail_n >= SC_SPP_FAIL_RESET) {
                int rnx = s_sc_ref_rtk->nx;
                memset(s_sc_ref_rtk->x, 0, rnx * sizeof(double));
                memset(s_sc_ref_rtk->P, 0, rnx * rnx * sizeof(double));
                s_ref_spp_fail_n = 0;
                sdr_log(3, "$LOG,%.3f,REF_PPP_KF_RESET (SPP failed %d epochs)",
                    time, SC_SPP_FAIL_RESET);
            }

            // Seed reference clock states when not yet initialized
            {
                int np = s_sc_ref_rtk->opt.dynamics ? 9 : 3;
                for (int i = 0; i < NSYS; i++) {
                    int ic = np + i;
                    if (s_sc_ref_rtk->x[ic] == 0.0 &&
                        s_sc_ref_rtk->sol.dtr[i] != 0.0) {
                        double dtr = i == 0 ? s_sc_ref_rtk->sol.dtr[0] :
                            s_sc_ref_rtk->sol.dtr[0] + s_sc_ref_rtk->sol.dtr[i];
                        s_sc_ref_rtk->x[ic] = CLIGHT * dtr;
                        s_sc_ref_rtk->P[ic + ic * s_sc_ref_rtk->nx] = 60.0 * 60.0;
                    }
                }
            }

            pppos(s_sc_ref_rtk, obs_ref, nobs_ref, pvt->nav);

            // Dynamics startup fix (same as main PPP)
            if (s_sc_ref_rtk->opt.dynamics) {
                double *rP = s_sc_ref_rtk->P;
                int rnx = s_sc_ref_rtk->nx;
                double rpstd = rnx > 0 ? sqrt(rP[0]+rP[1+rnx]+rP[2+2*rnx]) : 0.0;
                if (rpstd > 200.0) {
                    for (int j = 3; j < 9; j++) {
                        s_sc_ref_rtk->x[j] = 0.0;
                        for (int k = 0; k < rnx; k++) {
                            rP[j+k*rnx] = 0.0; rP[k+j*rnx] = 0.0;
                        }
                    }
                }
            }

            // Log reference PPP/SPP result (DOP from reference ssat — 4-unknown set)
            const sol_t *rsol = &s_sc_ref_rtk->sol;
            if (rsol->stat) {
                double ep[6], pos[3], rP9[9], Q[9], rdop[4] = {0};
                /* dtr[0] is in meters after pppos (PPP, stat=6), seconds after pntpos (SPP). */
                double rdtr = (rsol->stat == SOLQ_PPP)
                              ? rsol->dtr[0] / CLIGHT : rsol->dtr[0];
                time2epoch(timeadd(rsol->time, rdtr), ep);
                ecef2pos(rsol->rr, pos);
                rP9[0]=rsol->qr[0]; rP9[4]=rsol->qr[1]; rP9[8]=rsol->qr[2];
                rP9[1]=rP9[3]=rsol->qr[3];
                rP9[5]=rP9[7]=rsol->qr[4];
                rP9[2]=rP9[6]=rsol->qr[5];
                covenu(pos, rP9, Q);
                {
                    double razels[MAXSAT * 2];
                    int rns = 0;
                    for (int i = 0; i < MAXSAT; i++) {
                        if (s_sc_ref_rtk->ssat[i].azel[1] > 0.0) {
                            razels[rns*2]   = s_sc_ref_rtk->ssat[i].azel[0];
                            razels[rns*2+1] = s_sc_ref_rtk->ssat[i].azel[1];
                            rns++;
                        }
                    }
                    dops(rns, razels, 0.0, rdop);
                }
                sdr_log(3, "$REFPOS,%.3f,%.0f,%.0f,%.0f,%.0f,%.0f,%.3f,"
                    "%.9f,%.9f,%.3f,%d,%d,%.3f,%.3f,%.3f,%.9f,%.1f,%.1f,%.1f,%.1f,%.3e,"
                    "%.4f,%.4f,%.4f",
                    time, ep[0],ep[1],ep[2],ep[3],ep[4],ep[5],
                    pos[0]*R2D, pos[1]*R2D, pos[2], rsol->stat, rsol->ns,
                    SQRT(Q[4]), SQRT(Q[0]), SQRT(Q[8]), rdtr,
                    rdop[0], rdop[1], rdop[2], rdop[3], rsol->dtr[5],
                    rsol->rr[3], rsol->rr[4], rsol->rr[5]);
            }
        }

        // Prevent kinematic startup divergence.
        // In IFLC mode, PPP needs 4 dual-freq (L1+L2) satellites. L2CM takes
        // 30-60 s to lock on enough satellites. During that window, prnaccelh
        // process noise accumulates in the velocity state (no measurements to
        // correct it), causing the position to random-walk km away from truth.
        // By the time 4 L2CM sats lock, residuals >> maxinno and the KF can
        // never recover. Fix: zero vel/acc states each epoch while pos_std is
        // large (>200 m, i.e., PPP not yet converged), keeping the position
        // anchored to the SPP seed. Once PPP converges (pos_std falls below
        // 200 m), the dynamics model takes over normally.
        if (pvt->rtk->opt.dynamics) {
            double *P = pvt->rtk->P;
            int nx = pvt->rtk->nx;
            double pos_std = nx > 0 ? sqrt(P[0] + P[1+nx] + P[2+2*nx]) : 0.0;
            if (pos_std > 200.0) {
                for (int j = 3; j < 9; j++) {
                    pvt->rtk->x[j] = 0.0;
                    for (int k = 0; k < nx; k++) {
                        P[j + k*nx] = 0.0;
                        P[k + j*nx] = 0.0;
                    }
                }
            }
        }

        // Diagnostics: log KF position, clock, troposphere ZWD, LC ambiguity
        // convergence, and post-fit carrier residuals for PPP convergence monitoring.
        //
        // State vector layout (ppp.c macros, IONOOPT_IFLC, TROPOPT_ESTG):
        //   NP = dynamics ? 9 : 3   (pos/vel/acc)
        //   IC(0) = NP              (GPS clock, metres)
        //   IT     = NP + NSYS      (ZWD)
        //   NT     = 3              (ZWD + NS/EW gradients)
        //   NR     = NP + NSYS + NT (start of LC ambiguities)
        //   IB(s)  = NR + s - 1    (LC ambiguity for satellite s)
        {
            double *x = pvt->rtk->x;
            double *P = pvt->rtk->P;
            int nx = pvt->rtk->nx;
            int ndualfreq = 0;
            for (int i = 0; i < nobs; i++) {
                if (obs[i].code[1] && obs[i].P[1] != 0.0 && obs[i].L[1] != 0.0) ndualfreq++;
            }
            double pos_std = nx > 0 ? sqrt(P[0]+P[1+nx]+P[2+2*nx]) : -1;

            int np_i = pvt->rtk->opt.dynamics ? 9 : 3;
            double x_clk  = nx > np_i ? x[np_i] : 0.0;
            double clk_std = nx > np_i ? sqrt(P[np_i + np_i*nx]) : 0.0;

            /* troposphere ZWD: IT = NP + NSYS */
            int it = np_i + NSYS;
            int nt = pvt->rtk->opt.tropopt == TROPOPT_ESTG ? 3 :
                     pvt->rtk->opt.tropopt == TROPOPT_EST  ? 1 : 0;
            double zwd     = (nt > 0 && nx > it) ? x[it]              : 0.0;
            double zwd_std = (nt > 0 && nx > it) ? sqrt(P[it+it*nx])  : 0.0;

            /* LC ambiguities: IB(s,0) = NR + s - 1, NR = NP + NSYS + NT */
            int nr = np_i + NSYS + nt;
            double amb_var_sum = 0.0;
            double res_c_sq   = 0.0;
            int namb = 0, nconv = 0, nres = 0;
            for (int i = 0; i < MAXSAT; i++) {
                if (!pvt->rtk->ssat[i].vsat[0]) continue;
                int ib = nr + i;
                if (nx > ib && x[ib] != 0.0) {
                    amb_var_sum += P[ib + ib*nx];
                    namb++;
                    if (sqrt(P[ib + ib*nx]) < 0.5) nconv++;
                }
                double rc = pvt->rtk->ssat[i].resc[0];
                if (rc != 0.0) { res_c_sq += rc * rc; nres++; }
            }
            double amb_std = namb > 0 ? sqrt(amb_var_sum / namb) : 0.0;
            double res_c   = nres > 0 ? sqrt(res_c_sq   / nres)  : 0.0;

            sdr_log(3, "$LOG,%.3f,PPPOS_DBG"
                " kf_pos=%.1f,%.1f,%.1f pos_std=%.3f"
                " ndual=%d stat=%d clk=%.3f clk_std=%.3f"
                " zwd=%.4f zwd_std=%.4f amb_std=%.3f nconv=%d res_c=%.4f",
                time, x[0], x[1], x[2], pos_std, ndualfreq, pvt->rtk->sol.stat,
                x_clk, clk_std, zwd, zwd_std, amb_std, nconv, res_c);
        }

        // udbias_ppp() increments outc every epoch; update_stat() resets it but
        // only when stat==SOLQ_PPP. During convergence (stat=SOLQ_SINGLE), outc
        // keeps climbing and phase biases are wiped every maxout epochs even for
        // continuously tracked sats. Reset outc here for any sat with valid L2
        // obs so genuine outages (sat not in obs this epoch) still trigger resets.
        for (int i = 0; i < nobs; i++) {
            int valid_l1 = obs[i].code[0] && obs[i].P[0] != 0.0 && obs[i].L[0] != 0.0;
            int valid_l2 = obs[i].code[1] && obs[i].P[1] != 0.0 && obs[i].L[1] != 0.0;
            if ((sdr_ionoopt == IONOOPT_GRAPHIC) ? valid_l1 : valid_l2) {
                pvt->rtk->ssat[obs[i].sat - 1].outc[0] = 0;
            }
        }

        *pvt->sol = pvt->rtk->sol;
        memcpy(pvt->ssat, pvt->rtk->ssat, sizeof(ssat_t) * MAXSAT);

        if (sc_clk_ready) pvt->sol->dtr[5] = s_sc_clock_drift;
        if (pvt->sol->stat) {
            output_sol(pvt, time);
        } else if (spp_sol.stat) {
            // PPP has no solution yet (L2 not available or not converged).
            // Fall back to SPP so position output starts as soon as L1-only
            // observations are ready, rather than staying silent until L2 locks.
            *pvt->sol = spp_sol;
            if (sc_clk_ready) pvt->sol->dtr[5] = s_sc_clock_drift;
            output_sol(pvt, time);
        } else {
            /* Clear rr so the display shows no position rather than the
             * corrupted/seeded position from a failed solver. */
            pvt->sol->rr[0] = pvt->sol->rr[1] = pvt->sol->rr[2] = 0.0;
            update_azel(pvt->nav, pvt->sol, pvt->ssat);
            pvt->sol->ns = 0;
            sdr_log(3, "$LOG,%.3f,PPPOS NO SOLUTION", time);
        }
    } else {
        // SPP: single-point positioning with pseudorange
        prcopt_t opt = prcopt_default;
        opt.navsys |= SYS_GLO | SYS_GAL | SYS_QZS | SYS_CMP | SYS_IRN;
        opt.err[1] = opt.err[2] = STD_ERR;
        opt.ionoopt = sdr_ionoopt;
        opt.tropopt = TROPOPT_SAAS;
        opt.elmin   = sdr_el_mask * D2R;
        opt.maxgdop = sdr_maxgdop;
        opt.posopt[4] = 1; // RAIM-FDE
        if (sc_clk_ready) opt.clock_bias_fixed = 1; /* obs pre-corrected → 3-unknown */

        if (sdr_ps_sc_mode && !sc_clk_ready) {
            /* SC mode but no AOWR clock yet: suppress main solution */
            update_azel(pvt->nav, pvt->sol, pvt->ssat);
            pvt->sol->ns = 0;
            sdr_log(3, "$LOG,%.3f,PNTPOS SUPPRESSED (SC no clock)", time);
        } else {
            /* With clock_bias_fixed=1, pntpos fixes its internal clock state to
             * sol->dtr[0]*CLIGHT.  Must be set to clock_est (s) before the call;
             * without this x[3]=0 and residuals are off by ~clock_est*c (~3350 km). */
            if (sc_clk_ready) pvt->sol->dtr[0] = clock_est;

            if (pntpos(obs, nobs, pvt->nav, &opt, pvt->sol, NULL, pvt->ssat, msg)) {
                sdr_log(3, "$LOG,%.3f,SPP_SEED stat=%d dtr=%.9f msg=ok",
                    time, pvt->sol->stat, pvt->sol->dtr[0]);
                if (sc_clk_ready) pvt->sol->dtr[5] = s_sc_clock_drift;
                output_sol(pvt, time);
            } else {
                sdr_log(3, "$LOG,%.3f,SPP_SEED stat=0 dtr=0.000000000 msg=%s", time, msg);
                update_azel(pvt->nav, pvt->sol, pvt->ssat);
                pvt->sol->ns = 0;
            }
        }
    }
    pvt->nsat = pvt->obs->n;

    // SC reference solver: normal 4-unknown SPP on uncorrected obs for AOWR evaluation.
    // PPP mode uses its own reference pppos() above; this handles SPP mode only.
    if (run_ref_sol && sdr_pmode < PMODE_PPP_KINEMA) {
        prcopt_t ref_opt = prcopt_default;
        ref_opt.navsys |= SYS_GLO | SYS_GAL | SYS_QZS | SYS_CMP | SYS_IRN;
        ref_opt.err[1] = ref_opt.err[2] = STD_ERR;
        ref_opt.ionoopt = sdr_ionoopt;
        ref_opt.tropopt = TROPOPT_SAAS;
        ref_opt.elmin   = sdr_el_mask * D2R;
        ref_opt.maxgdop = sdr_maxgdop;
        ref_opt.posopt[4] = 1; /* RAIM-FDE */
        sol_t ref_sol = {0};
        double ref_azel[MAXSAT * 2] = {0}; /* azel from 4-unknown solve for DOP */
        int ref_stat = pntpos(obs_ref, nobs_ref, pvt->nav, &ref_opt, &ref_sol,
                              ref_azel, NULL, msg);
        sdr_log(3, "$LOG,%.3f,REF_SPP stat=%d dtr=%.9f msg=%s",
            time, ref_sol.stat, ref_sol.dtr[0], ref_stat ? "" : msg);
        if (ref_stat) {
            s_sc_ref_sol = ref_sol; /* keep for sdr_pvt_solstr() console display */
            double ep[6], pos[3], P[9], Q[9], rdop[4] = {0};
            double dtr_s = ref_sol.dtr[0];
            time2epoch(timeadd(ref_sol.time, dtr_s), ep);
            ecef2pos(ref_sol.rr, pos);
            P[0]=ref_sol.qr[0]; P[4]=ref_sol.qr[1]; P[8]=ref_sol.qr[2];
            P[1]=P[3]=ref_sol.qr[3]; P[5]=P[7]=ref_sol.qr[4]; P[2]=P[6]=ref_sol.qr[5];
            covenu(pos, P, Q);
            /* DOP from 4-unknown reference solve (independent of main clock-fixed DOP) */
            {
                double razels[MAXSAT * 2];
                int rns = 0;
                for (int i = 0; i < nobs_ref; i++) {
                    if (ref_azel[2*i+1] > ref_opt.elmin) {
                        razels[rns*2]   = ref_azel[2*i];
                        razels[rns*2+1] = ref_azel[2*i+1];
                        rns++;
                    }
                }
                dops(rns, razels, ref_opt.elmin, rdop);
            }
            sdr_log(3, "$REFPOS,%.3f,%.0f,%.0f,%.0f,%.0f,%.0f,%.3f,"
                "%.9f,%.9f,%.3f,%d,%d,%.3f,%.3f,%.3f,%.9f,%.1f,%.1f,%.1f,%.1f,%.3e,"
                "%.4f,%.4f,%.4f",
                time, ep[0],ep[1],ep[2],ep[3],ep[4],ep[5],
                pos[0]*R2D, pos[1]*R2D, pos[2], ref_sol.stat, ref_sol.ns,
                SQRT(Q[4]), SQRT(Q[0]), SQRT(Q[8]), dtr_s,
                rdop[0], rdop[1], rdop[2], rdop[3], ref_sol.dtr[5],
                ref_sol.rr[3], ref_sol.rr[4], ref_sol.rr[5]);
        }
    }

    // for debug
    double pos[3];
    ecef2pos(pvt->sol->rr, pos);
    trace(3, "%s %12.8f %13.8f %8.2f %d %2d/%2d DTR=%.1f %.1f %.1f %.1f (%s)\n",
        time_str(pvt->sol->time, 9), pos[0] * R2D, pos[1] * R2D, pos[2],
        pvt->sol->stat, pvt->sol->ns, pvt->nsat, pvt->sol->dtr[0] * 1e9,
        pvt->sol->dtr[1] * 1e9, pvt->sol->dtr[2] * 1e9, pvt->sol->dtr[3] * 1e9,
        msg);
    for (int i = 0; i < MAXSAT; i++) {
        ssat_t *ssat = pvt->ssat + i;
        if (ssat->azel[1] <= 0.0) continue;
        char sat[16];
        satno2id(i+1, sat);
        trace(3, "%s %d %4.1f %5.1f %4.1f %12.3f\n", sat, ssat->vs,
            ssat->snr[0] * SNR_UNIT, ssat->azel[0] * R2D, ssat->azel[1] * R2D,
            ssat->resp[0]);
    }
    s_rx_clock_s  = (pvt->sol->stat == SOLQ_PPP) ?
        pvt->sol->dtr[0] / CLIGHT : pvt->sol->dtr[0];
    s_rx_sol_valid = (pvt->sol->stat > 0);
    // Write GS ring-buffer at every PVT epoch once PS has been observed,
    // matching gnss-sdr write_clock_difference() which persists after PS loss.
    if (sdr_ps_gs_mode && s_rx_sol_valid && s_gs_ps_ever_obs) {
        int week;
        double tow = time2gpst(pvt->sol->time, &week);
        write_aowr_gs(tow - s_gs_last_dt_aowr, -s_gs_last_dt_aowr + s_rx_clock_s);
    }
}

// resolve msec ambiguity in pseudorange ---------------------------------------
static void res_obs_amb(obs_t *obs, int sys, uint8_t code, double sec)
{
    for (int i = 0; i < obs->n; i++) {
        obsd_t *data = obs->data + i;
        if (!(satsys(data->sat, NULL) & sys)) continue;
        
        for (int j = 0; j < NFREQ + NEXOBS; j++) {
            if (data->code[j] != code) continue;
            int k;
            for (k = 0; k < NFREQ + NEXOBS; k++) {
                if (!data->code[k] || data->code[k] == code ||
                    data->code[k] == CODE_L5Q || data->code[k] == CODE_L5P) {
                    continue;
                }
                double tau1 = data->P[j] / CLIGHT, tau2 = data->P[k] / CLIGHT;
                double tau3 = floor(tau2 / sec) * sec + fmod(tau1, sec);
                if      (tau3 < tau2 - sec / 2.0) tau3 += sec;
                else if (tau3 > tau2 + sec / 2.0) tau3 -= sec;
                data->P[j] = CLIGHT * tau3;
                break;
            }
            if (k >= NFREQ + NEXOBS) {
                data->P[j] = 0.0; // set invalid if unresolved
            }
        }
    }
}

//------------------------------------------------------------------------------
//  Update PVT solution.
//
//  args:
//      pvt      (IO) SDR PVT
//      ix       (I)  received IF data cycle (cyc)
//
//  returns:
//      none
//
void sdr_pvt_udsol(sdr_pvt_t *pvt, int64_t ix)
{
    sdr_mutex_lock(&pvt->mtx);
    
    if (pvt->ix > 0 && (pvt->nch >= pvt->rcv->nch ||
        ix >= pvt->ix + (int)(sdr_lag_epoch / SDR_CYC))) {
        
        // sort obs data (ordering for the merge step in update_sol)
        sortobs(pvt->obs);
        if (pvt->obs->n > 0) pvt->count[1]++;
        
        // update PVT solution
        update_sol(pvt);
        
        sdr_rcv_array_calib(pvt->rcv, pvt->obs->data, pvt->obs->n,
            pvt->nav, pvt->sol->rr);
        
        // solution latency (s)
        pvt->latency = (ix - pvt->ix) * SDR_CYC;
        
        // set next epoch time and cycle
        pvt->time = timeadd(pvt->time, sdr_epoch);
        pvt->ix += (int)(sdr_epoch / SDR_CYC);
        pvt->nch = pvt->obs->n = 0; 
        
        // adjust epoch cycle within 20 ms
        // pntpos (stat=SINGLE) stores dtr in seconds; pppos (stat=PPP) stores in meters.
        // Same rule as out_log_pos(): check sol->stat, not sdr_pmode.
        // Skip for GS side: epoch adjustment would corrupt AOWR timing fed to SC.
        // Skip for SC side once AOWR is providing the clock: pvt->ix must not shift
        // because the AOWR clock manages the receiver epoch externally.
        if (pvt->sol->stat && !sdr_ps_gs_mode &&
                !(sdr_ps_sc_mode && s_sc_hist_n > 0)) {
            double dtr_s = pvt->sol->dtr[0];
            if (pvt->sol->stat == SOLQ_PPP)
                dtr_s /= CLIGHT;
            double dtr = ROUND(dtr_s / 0.02) * 0.02;
            if (fabs(dtr) > 0.01) {
                pvt->ix += (int)(dtr / SDR_CYC);
                sdr_log(3, "$LOG,%.3f,PVT EPOCH ADJUSTED (DT=%.3fs)",
                    pvt->ix * SDR_CYC, dtr);
                /* Epoch shift changes the local time reference by dtr seconds.
                 * All AOWR clock offsets (which measure receiver_local - GPS_time)
                 * must be reduced by the same amount to remain coherent. */
                if (sdr_ps_sc_mode) {
                    s_sc_dt_i         -= dtr;
                    s_sc_last_dt_aowr -= dtr;
                    for (int k = 0; k < s_sc_hist_n; k++) {
                        s_sc_clk_hist[k] -= dtr;
                        /* time2gpst(obs[0].time) advances by dtr after each ix
                         * correction; tow_hist must track the same shift so dt
                         * (current_tow - tow_hist[0]) stays smooth. */
                        s_sc_tow_hist[k] += dtr;
                    }
                }
            }
        }
    }
    sdr_mutex_unlock(&pvt->mtx);
}

//------------------------------------------------------------------------------
//  Get PVT solution string.
//
//  args:
//      pvt      (I)  SDR PVT
//      buff     (IO) PVT solution string buffer
//      size     (I)  size of string buffer
//
//  returns:
//      none
//
void sdr_pvt_solstr(sdr_pvt_t *pvt, char *buff, int size)
{
    static const char *solq[] = {"---","FIX","FLT","SBS","DGP","SPP","PPP","DR"};
    char tstr[32] = "", nstr[16] = "";
    double pos[3] = {0};
    int stat = 0;

    sdr_mutex_lock(&pvt->mtx);

    if (norm(pvt->sol->rr, 3) > 1e-6) {
        time2str(pvt->sol->time, tstr, 1);
        ecef2pos(pvt->sol->rr, pos);
        stat = pvt->sol->stat;
    } else {
        time2str(pvt->time, tstr, 1);
    }
    sdr_mutex_unlock(&pvt->mtx);

    tstr[4] = tstr[7] = '-';
    snprintf(nstr, sizeof(nstr), "%d/%d", pvt->sol->ns, pvt->nsat);

    if (sdr_ps_sc_mode) {
        /* ps_sc mode: show AOWR-corrected (SC) and uncorrected (REF) solutions */
        char rtstr[32] = "", rnstr[16] = "";
        double rpos[3] = {0};
        int rstat = 0;
        const sol_t *rsol = s_sc_ref_rtk ? &s_sc_ref_rtk->sol : &s_sc_ref_sol;
        if (norm(rsol->rr, 3) > 1e-6) {
            time2str(rsol->time, rtstr, 1);
            ecef2pos(rsol->rr, rpos);
            rstat = rsol->stat;
        } else {
            time2str(pvt->time, rtstr, 1);
        }
        rtstr[4] = rtstr[7] = '-';
        snprintf(rnstr, sizeof(rnstr), "%d/%d", rsol->ns, pvt->nsat);
        snprintf(buff, size,
            "SC %21s %12.8f %13.8f %9.3f %5s %s\n"
            "REF%21s %12.8f %13.8f %9.3f %5s %s",
            tstr,  pos[0]*R2D,  pos[1]*R2D,  pos[2],  nstr,
            solq[stat  >= 0 && stat  < 8 ? stat  : 0],
            rtstr, rpos[0]*R2D, rpos[1]*R2D, rpos[2], rnstr,
            solq[rstat >= 0 && rstat < 8 ? rstat : 0]);
    } else {
        snprintf(buff, size, "%21s %12.8f %13.8f %9.3f %5s %s", tstr, pos[0] * R2D,
            pos[1] * R2D, pos[2], nstr,
            solq[stat >= 0 && stat < 8 ? stat : 0]);
    }
}
