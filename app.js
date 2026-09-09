const sigma = 5.670374419e-8;

const defaults = {
  solarConstant: 1361,
  albedo: 0.30,
  emissivity: 0.80
};

const solarSlider = document.getElementById("solarSlider");
const solarInput = document.getElementById("solarInput");
const albedoSlider = document.getElementById("albedoSlider");
const albedoInput = document.getElementById("albedoInput");
const emissivitySlider = document.getElementById("emissivitySlider");
const emissivityInput = document.getElementById("emissivityInput");

function linkPair(slider, input, decimals) {
  slider.addEventListener("input", () => {
    input.value = Number(slider.value).toFixed(decimals);
  });

  input.addEventListener("input", () => {
    const value = Number(input.value);
    if (Number.isFinite(value)) {
      slider.value = value;
    }
  });
}

linkPair(solarSlider, solarInput, 0);
linkPair(albedoSlider, albedoInput, 2);
linkPair(emissivitySlider, emissivityInput, 2);

function runRadiationModel(S0, alpha, epsilon) {
  const SW_IN = S0 / 4;
  const SWrefl = alpha * SW_IN;
  const SWabs = (1 - alpha) * SW_IN;

  const denom = sigma * (1 - 0.5 * epsilon);
  const Tsurf_K = Math.pow(SWabs / denom, 0.25);
  const Tsurf_C = Tsurf_K - 273.1;

  let Tatm_K = NaN;
  if (epsilon !== 0) {
    Tatm_K = Math.pow(Math.pow(Tsurf_K, 4) / 2, 0.25);
  }
  const Tatm_C = Tatm_K - 273.1;

  const LWout_surf = sigma * Math.pow(Tsurf_K, 4);
  const LWsurf_atm = epsilon * LWout_surf;
  const LWsurf_space = (1 - epsilon) * LWout_surf;

  let LWatm_surf = 0;
  let LWatm_space = 0;

  if (epsilon !== 0) {
    LWatm_surf = epsilon * sigma * Math.pow(Tatm_K, 4);
    LWatm_space = epsilon * sigma * Math.pow(Tatm_K, 4);
  }

  const Enet_surf = SWabs + LWatm_surf - LWout_surf;
  const Enet_atm = LWsurf_atm - LWatm_surf - LWatm_space;
  const Enet_space = SW_IN - SWrefl - LWsurf_space - LWatm_space;

  return {
    incomingSW: SW_IN,
    reflectedSW: SWrefl,
    absorbedSW: SWabs,
    surfaceLW: LWout_surf,
    atmosAbsorbedLW: LWsurf_atm,
    lwWindow: LWsurf_space,
    atmosDownwardLW: LWatm_surf,
    atmosUpwardLW: LWatm_space,
    surfaceTemperatureC: Tsurf_C,
    airTemperatureC: Tatm_C,
    surfaceBalance: Enet_surf,
    atmosphereBalance: Enet_atm,
    toaBalance: Enet_space
  };
}

function validateInputs(S0, alpha, epsilon) {
  if (!Number.isFinite(S0) || S0 < 800 || S0 > 1600) {
    throw new Error("Solar constant must be between 800 and 1600 W m⁻².");
  }
  if (!Number.isFinite(alpha) || alpha < 0 || alpha > 1) {
    throw new Error("Earth albedo must be between 0 and 1.");
  }
  if (!Number.isFinite(epsilon) || epsilon < 0 || epsilon > 1) {
    throw new Error("Atmospheric LW emissivity must be between 0 and 1.");
  }
}

function formatFlux(value) {
  return Math.round(value).toString();
}

function calculateAndDisplay() {
  const status = document.getElementById("status");

  try {
    const S0 = Number(solarInput.value);
    const alpha = Number(albedoInput.value);
    const epsilon = Number(emissivityInput.value);

    validateInputs(S0, alpha, epsilon);
    const r = runRadiationModel(S0, alpha, epsilon);

    document.getElementById("incomingSW").textContent = formatFlux(r.incomingSW);
    document.getElementById("reflectedSW").textContent = formatFlux(r.reflectedSW);
    document.getElementById("absorbedSW").textContent = formatFlux(r.absorbedSW);

    document.getElementById("surfaceLW").textContent = formatFlux(r.surfaceLW);
    document.getElementById("atmosAbsorbedLW").textContent = formatFlux(r.atmosAbsorbedLW);
    document.getElementById("lwWindow").textContent = formatFlux(r.lwWindow);
    document.getElementById("atmosDownwardLW").textContent = formatFlux(r.atmosDownwardLW);
    document.getElementById("atmosUpwardLW").textContent = formatFlux(r.atmosUpwardLW);

    document.getElementById("surfaceTemp").textContent =
      `Ts: ${r.surfaceTemperatureC.toFixed(1)} °C`;

    document.getElementById("airTemp").textContent =
      Number.isNaN(r.airTemperatureC)
        ? "Tatm: N/A"
        : `Ta: ${r.airTemperatureC.toFixed(1)} °C`;

    status.textContent =
      `Calculated successfully\nTOA imbalance: ${r.toaBalance.toExponential(3)} W m⁻²`;
    status.className = "status ok";
  } catch (err) {
    status.textContent = err.message;
    status.className = "status error";
  }
}

function resetInputs() {
  solarSlider.value = defaults.solarConstant;
  solarInput.value = defaults.solarConstant.toFixed(0);

  albedoSlider.value = defaults.albedo;
  albedoInput.value = defaults.albedo.toFixed(2);

  emissivitySlider.value = defaults.emissivity;
  emissivityInput.value = defaults.emissivity.toFixed(2);

  calculateAndDisplay();
}

document.getElementById("calculateButton").addEventListener("click", calculateAndDisplay);
document.getElementById("resetButton").addEventListener("click", resetInputs);

resetInputs();
