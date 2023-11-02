from ixdat import Measurement


ec = Measurement.read(
    "../data/01/iridium_butterfly_short_CVA.mpt", reader="biologic"
)

ec.calibrate(RE_vs_RHE=0.720, A_el=0.196)
cv = ec.as_cv()

cv.redefine_cycle(start_potential=0.3, redox=True)

ax = cv.plot_cycles()
fig = ax.get_figure()
fig.set_figwidth(6)
fig.set_figheight(4.5)
fig.tight_layout()
fig.savefig("demo_01.png", dpi=600)
