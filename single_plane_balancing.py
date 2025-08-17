import streamlit as st
import math, cmath
import numpy as np
import matplotlib.pyplot as plt
from io import BytesIO
from reportlab.pdfgen import canvas
from reportlab.lib.units import cm

# ----------------------
# Helpers
# ----------------------
def norm_deg(d):
    return d % 360

def polar_to_rect(mag, angle_deg):
    angle_rad = math.radians(angle_deg)
    return complex(mag * math.cos(angle_rad), mag * math.sin(angle_rad))

def rect_to_polar(z: complex):
    mag = abs(z)
    ang = (math.degrees(cmath.phase(z)) + 360.0) % 360.0
    return mag, ang

def complex_from_mass(mag_g, angle_deg):
    return mag_g * np.exp(1j * np.deg2rad(angle_deg))

# ----------------------
# PDF builder
# ----------------------
def build_pdf_report(title, history, baseline=None, trial=None):
    buf = BytesIO()
    c = canvas.Canvas(buf)
    c.setFont("Helvetica-Bold", 14)
    c.drawString(2*cm, 28*cm, f"Rotor Balancing Report: {title}")
    y = 26*cm

    if baseline:
        c.setFont("Helvetica-Bold", 12)
        c.drawString(2*cm, y, "Baseline Vibration:")
        y -= 0.8*cm
        for plane, val in baseline.items():
            c.setFont("Helvetica", 12)
            c.drawString(3*cm, y, f"{plane}: {val['mag']:.2f} @ {val['ang']:.1f}°")
            y -= 0.6*cm
        y -= 0.4*cm

    if trial:
        c.setFont("Helvetica-Bold", 12)
        c.drawString(2*cm, y, "Trial Masses:")
        y -= 0.8*cm
        for plane, val in trial.items():
            c.setFont("Helvetica", 12)
            c.drawString(3*cm, y, f"{plane}: {val['mass']:.2f} g @ {val['ang']:.1f}°")
            y -= 0.6*cm
        y -= 0.4*cm

    for step in history:
        c.setFont("Helvetica", 12)
        desc = step.get("desc","")
        c.drawString(2*cm, y, f"Step {step['step']+1}: {desc}")
        if "details" in step:
            for key, val in step["details"].items():
                y -= 0.6*cm
                c.drawString(3*cm, y, f"{key}: {val}")
        y -= 1.2*cm
        if y < 2*cm:
            c.showPage()
            y = 28*cm

    c.showPage()
    c.save()
    buf.seek(0)
    return buf

# ----------------------
# Reset function
# ----------------------
def reset_app():
    for key in list(st.session_state.keys()):
        del st.session_state[key]
    st.rerun()

# ----------------------
# Single-plane flow
# ----------------------
def single_plane_flow():
    st.header("⚙️ Single Plane Balancing")
    if "sp_step" not in st.session_state: st.session_state.sp_step = 0
    if "sp" not in st.session_state: st.session_state.sp = {}
    if "sp_history" not in st.session_state: st.session_state.sp_history = []
    if "sp_baseline" not in st.session_state: st.session_state.sp_baseline = None
    if "sp_trial" not in st.session_state: st.session_state.sp_trial = None

    if st.session_state.sp_step == 0:
        with st.form("sp_trial_form"):
            st.subheader("Step 1: Enter Trial Data")
            target = st.number_input("Target amplitude", value=0.5, step=0.1)
            V0_mag = st.number_input("Baseline amplitude", value=4.0, step=0.1)
            V0_ang = norm_deg(st.number_input("Baseline phase (deg)", value=20.0, step=1.0))
            trial_mass = st.number_input("Trial mass (g)", value=10.0, step=0.1)
            trial_angle = norm_deg(st.number_input("Trial mass angle (deg)", value=0.0, step=1.0))
            V1_mag = st.number_input("Amplitude with trial mass", value=1.5, step=0.1)
            V1_ang = norm_deg(st.number_input("Phase with trial mass (deg)", value=-160.0, step=1.0))
            submit = st.form_submit_button("➡️ Submit & Continue")
        if submit:
            V0 = polar_to_rect(V0_mag, V0_ang)
            V1 = polar_to_rect(V1_mag, V1_ang)
            m_trial = complex_from_mass(trial_mass, trial_angle)
            a = (V1 - V0)/m_trial
            st.session_state.sp = {"target": target, "trial_angle": trial_angle, "a": a, "current_vib": V0}
            st.session_state.sp_step = 1
            st.session_state.sp_baseline = {"Baseline": {"mag": V0_mag, "ang": V0_ang}}
            st.session_state.sp_trial = {"Trial mass": {"mass": trial_mass, "ang": trial_angle}}
            st.rerun()
    else:
        data = st.session_state.sp
        current_vib = data["current_vib"]
        a = data["a"]
        trial_angle = data["trial_angle"]
        m_corr = -current_vib / a
        m_corr_mag, m_corr_ang_rel = rect_to_polar(m_corr)
        final_angle = norm_deg(trial_angle + m_corr_ang_rel)
        st.subheader(f"Step {st.session_state.sp_step + 1}: Correction")
        st.success(f"Suggested correction: **{m_corr_mag:.3f} g** at **{final_angle:.1f}°**")

        fig, ax = plt.subplots(figsize=(5,5))
        ax.quiver(0,0,current_vib.real,current_vib.imag,angles='xy',scale_units='xy',scale=1,label="Current Vib")
        ax.quiver(0,0,m_corr.real,m_corr.imag,angles='xy',scale_units='xy',scale=1,label="Correction Mass")
        max_val = max(abs(current_vib), m_corr_mag)*1.4 if max(abs(current_vib), m_corr_mag)>0 else 1
        ax.set_xlim(-max_val, max_val); ax.set_ylim(-max_val, max_val)
        ax.axhline(0,color='black',linewidth=0.5); ax.axvline(0,color='black',linewidth=0.5)
        ax.set_aspect('equal','box'); ax.grid(True); ax.legend(); ax.set_title("Single-plane Phasor")
        st.pyplot(fig)

        exact = st.radio("Did you apply the exact suggested correction?", ["Yes","No"], key=f"sp_exact_{st.session_state.sp_step}")
        with st.form(f"sp_iter_form_{st.session_state.sp_step}"):
            if exact=="No":
                applied_mass = st.number_input("Applied mass (g)", value=0.0, step=0.1, key=f"sp_mass_{st.session_state.sp_step}")
                applied_angle = norm_deg(st.number_input("Applied angle (deg)", value=0.0, step=1.0, key=f"sp_ang_{st.session_state.sp_step}"))
            else:
                applied_mass, applied_angle = m_corr_mag, final_angle
            new_mag = st.number_input("New amplitude after correction", value=0.0, step=0.1, key=f"sp_newmag_{st.session_state.sp_step}")
            new_ang = norm_deg(st.number_input("New phase after correction (deg)", value=0.0, step=1.0, key=f"sp_newang_{st.session_state.sp_step}"))
            go = st.form_submit_button("🔄 Next")
        if go:
            st.session_state.sp_history.append({
                "step": st.session_state.sp_step,
                "desc": f"Applied correction",
                "details": {
                    "Applied mass": f"{applied_mass:.2f} g @ {applied_angle:.1f}°",
                    "Resulting vibration": f"{new_mag:.2f} @ {new_ang:.1f}°"
                }
            })
            st.session_state.sp["current_vib"] = polar_to_rect(new_mag,new_ang)
            if abs(st.session_state.sp["current_vib"]) <= data["target"]:
                st.success("✅ Vibration below target. Balancing complete!")
                pdf_buf = build_pdf_report("Single Plane", st.session_state.sp_history,
                                           baseline=st.session_state.sp_baseline,
                                           trial=st.session_state.sp_trial)
                col1, col2 = st.columns([1,1])
                with col1: st.download_button("📄 Download PDF", pdf_buf, file_name="SinglePlane_Balancing.pdf", mime="application/pdf")
                with col2: st.button("🔄 Reset", on_click=reset_app)
            else:
                st.session_state.sp_step +=1
                st.rerun()

# ----------------------
# Two-plane flow (integrated)
# ----------------------
def two_plane_flow():
    st.header("⚙️ Two Plane Balancing — Simultaneous")
    if "tp_step" not in st.session_state: st.session_state.tp_step = 0
    if "tp" not in st.session_state: st.session_state.tp = {}
    if "tp_history" not in st.session_state: st.session_state.tp_history = []
    if "tp_baseline" not in st.session_state: st.session_state.tp_baseline = None
    if "tp_trial" not in st.session_state: st.session_state.tp_trial = None

    if st.session_state.tp_step == 0:
        with st.form("tp_trial_form"):
            st.subheader("Step 1: Enter Baseline & Trial Data")
            target = st.number_input("Target amplitude for both planes", value=0.5, step=0.1)

            st.markdown("**Baseline vibrations**")
            V0A_mag = st.number_input("Plane A baseline amplitude", value=4.0, step=0.1)
            V0A_ang = norm_deg(st.number_input("Plane A baseline phase (deg)", value=20.0, step=1.0))
            V0B_mag = st.number_input("Plane B baseline amplitude", value=3.0, step=0.1)
            V0B_ang = norm_deg(st.number_input("Plane B baseline phase (deg)", value=50.0, step=1.0))

            st.markdown("**Trial on Plane A**")
            tA_mass = st.number_input("Trial mass @ Plane A (g)", value=10.0, step=0.1)
            tA_angle = norm_deg(st.number_input("Trial angle @ Plane A (deg)", value=0.0, step=1.0))
            V1A_mag = st.number_input("Plane A amplitude with trial A", value=2.0, step=0.1)
            V1A_ang = norm_deg(st.number_input("Plane A phase with trial A (deg)", value=160.0, step=1.0))
            V1B_mag = st.number_input("Plane B amplitude with trial A", value=2.5, step=0.1)
            V1B_ang = norm_deg(st.number_input("Plane B phase with trial A (deg)", value=120.0, step=1.0))

            st.markdown("**Trial on Plane B**")
            tB_mass = st.number_input("Trial mass @ Plane B (g)", value=10.0, step=0.1)
            tB_angle = norm_deg(st.number_input("Trial angle @ Plane B (deg)", value=90.0, step=1.0))
            V2A_mag = st.number_input("Plane A amplitude with trial B", value=1.8, step=0.1)
            V2A_ang = norm_deg(st.number_input("Plane A phase with trial B (deg)", value=140.0, step=1.0))
            V2B_mag = st.number_input("Plane B amplitude with trial B", value=1.5, step=0.1)
            V2B_ang = norm_deg(st.number_input("Plane B phase with trial B (deg)", value=100.0, step=1.0))

            submit = st.form_submit_button("➡️ Submit & Compute Influence Matrix")
        if submit:
            V0A = polar_to_rect(V0A_mag,V0A_ang)
            V0B = polar_to_rect(V0B_mag,V0B_ang)
            V1A = polar_to_rect(V1A_mag,V1A_ang)
            V1B = polar_to_rect(V1B_mag,V1B_ang)
            V2A = polar_to_rect(V2A_mag,V2A_ang)
            V2B = polar_to_rect(V2B_mag,V2B_ang)
            mA = complex_from_mass(tA_mass,tA_angle)
            mB = complex_from_mass(tB_mass,tB_angle)

            A = np.array([[ (V1A-V0A)/mA , (V2A-V0A)/mB ],
                          [ (V1B-V0B)/mA , (V2B-V0B)/mB ]], dtype=complex)

            st.session_state.tp = {"target": target, "trial_angles":(tA_angle,tB_angle), "A":A, "current":np.array([V0A,V0B],dtype=complex)}
            st.session_state.tp_step = 1
            st.session_state.tp_baseline = {"Plane A": {"mag": V0A_mag, "ang": V0A_ang},
                                            "Plane B": {"mag": V0B_mag, "ang": V0B_ang}}
            st.session_state.tp_trial = {"Trial A": {"mass": tA_mass, "ang": tA_angle},
                                         "Trial B": {"mass": tB_mass, "ang": tB_angle}}
            st.rerun()

    else:
        tp = st.session_state.tp
        A = tp["A"]
        V = tp["current"]
        tA_angle, tB_angle = tp["trial_angles"]
        target = tp["target"]

        try:
            Ainv = np.linalg.inv(A)
        except np.linalg.LinAlgError:
            st.error("❌ Influence matrix not invertible")
            return

        m_vec = -Ainv @ V
        mA_c, mB_c = m_vec[0], m_vec[1]
        mA_mag, mA_ang_rel = rect_to_polar(mA_c)
        mB_mag, mB_ang_rel = rect_to_polar(mB_c)
        mA_final_angle = norm_deg(tA_angle + mA_ang_rel)
        mB_final_angle = norm_deg(tB_angle + mB_ang_rel)

        st.subheader(f"Step {st.session_state.tp_step + 1}: Simultaneous Correction")
        st.success(f"Plane A: {mA_mag:.3f} g @ {mA_final_angle:.1f}°")
        st.success(f"Plane B: {mB_mag:.3f} g @ {mB_final_angle:.1f}°")

        exact = st.radio("Did you apply exact corrections?", ["Yes","No"], key=f"tp_exact_{st.session_state.tp_step}")
        with st.form(f"tp_iter_form_{st.session_state.tp_step}"):
            if exact=="No":
                a_mass = st.number_input("Applied mass @ Plane A (g)", value=0.0, step=0.1, key=f"tp_A_mass_{st.session_state.tp_step}")
                a_angle = norm_deg(st.number_input("Applied angle @ Plane A (deg)", value=0.0, step=1.0, key=f"tp_A_ang_{st.session_state.tp_step}"))
                b_mass = st.number_input("Applied mass @ Plane B (g)", value=0.0, step=0.1, key=f"tp_B_mass_{st.session_state.tp_step}")
                b_angle = norm_deg(st.number_input("Applied angle @ Plane B (deg)", value=0.0, step=1.0, key=f"tp_B_ang_{st.session_state.tp_step}"))
            else:
                a_mass, a_angle = mA_mag, mA_final_angle
                b_mass, b_angle = mB_mag, mB_final_angle

            newA_amp = st.number_input("New amplitude @ Plane A", value=0.0, step=0.1, key=f"tp_newA_amp_{st.session_state.tp_step}")
            newA_ang = norm_deg(st.number_input("New phase @ Plane A", value=0.0, step=1.0, key=f"tp_newA_ang_{st.session_state.tp_step}"))
            newB_amp = st.number_input("New amplitude @ Plane B", value=0.0, step=0.1, key=f"tp_newB_amp_{st.session_state.tp_step}")
            newB_ang = norm_deg(st.number_input("New phase @ Plane B", value=0.0, step=1.0, key=f"tp_newB_ang_{st.session_state.tp_step}"))
            go = st.form_submit_button("🔄 Next")
        if go:
            st.session_state.tp_history.append({"step":st.session_state.tp_step,
                                               "desc": f"Applied correction",
                                               "details": {
                                                   "Plane A": f"{a_mass:.2f} g @ {a_angle:.1f}° → {newA_amp:.2f} @ {newA_ang:.1f}°",
                                                   "Plane B": f"{b_mass:.2f} g @ {b_angle:.1f}° → {newB_amp:.2f} @ {newB_ang:.1f}°"
                                               }})
            st.session_state.tp["current"] = np.array([polar_to_rect(newA_amp,newA_ang),
                                                       polar_to_rect(newB_amp,newB_ang)],dtype=complex)
            if abs(st.session_state.tp["current"][0]) <= target and abs(st.session_state.tp["current"][1]) <= target:
                st.success("✅ Two-plane balancing complete!")
                pdf_buf = build_pdf_report("Two Plane", st.session_state.tp_history,
                                           baseline=st.session_state.tp_baseline,
                                           trial=st.session_state.tp_trial)
                col1,col2=st.columns([1,1])
                with col1: st.download_button("📄 Download PDF", pdf_buf, file_name="TwoPlane_Balancing.pdf", mime="application/pdf")
                with col2: st.button("🔄 Reset", on_click=reset_app)
            else:
                st.session_state.tp_step += 1
                st.rerun()

# ----------------------
# Main App
# ----------------------
st.set_page_config(page_title="Rotor Balancing", layout="centered")
st.title("🌀 Rotor Balancing Calculator")

if "mode" not in st.session_state: st.session_state.mode="Single Plane"
if "started" not in st.session_state: st.session_state.started=False

mode_choice = st.radio("Select Balancing Mode:", ["Single Plane","Two Plane"],
                       index=0 if st.session_state.mode=="Single Plane" else 1, horizontal=True)
st.session_state.mode = mode_choice

if not st.session_state.started:
    with st.form("start_form"):
        st.subheader("Welcome! Choose balancing type and start")
        start_btn = st.form_submit_button("🚀 Start")
    if start_btn:
        st.session_state.started=True
        st.session_state.pop("sp_step",None); st.session_state.pop("tp_step",None)
        st.rerun()
else:
    st.caption(f"Mode: **{st.session_state.mode}**")
    if st.session_state.mode=="Single Plane":
        single_plane_flow()
    else:
        two_plane_flow()
