NO_TIME_CUT = 9999.0
TIME_CUT_STEP = 0.01
MAX_TIME_CUT = 0.50
SIGNAL_EFFICIENCIES = (0.90, 0.95)
PLOT_HLT_DEFAULT = False


def timing_cuts():
  count = round(MAX_TIME_CUT / TIME_CUT_STEP)
  return [round(TIME_CUT_STEP * (index + 1), 3) for index in range(count)] + [NO_TIME_CUT]
