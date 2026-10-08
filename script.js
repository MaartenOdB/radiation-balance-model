"use strict";

/*
===============================================================
RADIATION BALANCE MODEL
JavaScript version of radiationBalanceModel.m
===============================================================
*/


// =============================================================
// 1. GREENHOUSE-GAS FORCING
//    Meinshausen et al. (2020) formulation
// =============================================================

function forcingGHGMeinshausen(C, M, N) {

    // Preindustrial concentrations
    const Cpi = 278.3;     // CO2 [ppm]
    const Mpi = 729.2;     // CH4 [ppb]
    const Npi = 270.1;     // N2O [ppb]

    // AR6 coefficients

    // CO2
    const a1 = -2.48e-7;
    const b1 =  7.59e-4;
    const c1 = -2.15e-3;
    const d1 =  5.2488;
    const e1 =  1.05;

    // CH4
    const a2 = -8.96e-5;
    const b2 = -1.25e-4;
    const d2 =  0.045194;
    const e2 =  0.86;

    // N2O
    const a3 = -3.42e-4;
    const b3 =  2.546e-4;
    const c3 = -2.44e-4;
    const d3 =  0.12173;
    const e3 =  1.07;


    // ---------------------------------------------------------
    // CO2 forcing
    // ---------------------------------------------------------

    const C_alpha_max = Cpi - b1 / (2 * a1);

    let alpha_CO2;

    if (C <= Cpi) {

        alpha_CO2 = d1;

    } else if (C <= C_alpha_max) {

        alpha_CO2 =
            d1 +
            a1 * Math.pow(C - Cpi, 2) +
            b1 * (C - Cpi);

    } else {

        alpha_CO2 =
            d1 -
            Math.pow(b1, 2) / (4 * a1);
    }

    const F_CO2 =
        (alpha_CO2 + c1 * Math.sqrt(N)) *
        Math.log(C / Cpi) *
        e1;


    // ---------------------------------------------------------
    // CH4 forcing
    // ---------------------------------------------------------

    const F_CH4 =
        (a2 * Math.sqrt(M) +
         b2 * Math.sqrt(N) +
         d2) *
        (Math.sqrt(M) - Math.sqrt(Mpi)) *
        e2;


    // ---------------------------------------------------------
    // N2O forcing
    // ---------------------------------------------------------

    const F_N2O =
        (a3 * Math.sqrt(C) +
         b3 * Math.sqrt(N) +
         c3 * Math.sqrt(M) +
         d3) *
        (Math.sqrt(N) - Math.sqrt(Npi)) *
        e3;


    return F_CO2 + F_CH4 + F_N2O;
}



// =============================================================
// 2. RADIATION BALANCE MODEL
// =============================================================

function radiationBalanceModel(
    solarConstant,
    albedoAtmosphere,
    albedoSurface,
    CO2,
    CH4,
    N2O
) {

    // =========================================================
    // CONSTANTS AND REFERENCE VALUES
    // =========================================================

    const sigma = 5.670374419e-8;

    // Fixed atmospheric SW absorptivity
    const SWabsAtmosphere = 0.22;

    // Reference climate
    const S0  = 1361;
    const aa0 = 0.23;
    const as0 = 0.12;

    const CO2_0 = 424;
    const CH4_0 = 1935;
    const N2O_0 = 339;

    const T0_C = 15.0;
    const T0_K = T0_C + 273.15;

    // Reference turbulent fluxes
    const H0  = 21;
    const LE0 = 82;

    // Temperature dependence of turbulent fluxes
    const gammaH  = 0.01;
    const gammaLE = 0.025;

    // Effective net climate feedback parameter
    const lambda = 1.30;

    // Reference atmospheric window
    const LWwindow0 = 40;


    // =========================================================
    // CHECK INPUTS
    // =========================================================

    if (solarConstant < 0) {
        throw new Error("Solar constant cannot be negative.");
    }

    if (albedoAtmosphere < 0 || albedoAtmosphere > 1) {
        throw new Error(
            "Atmospheric albedo must be between 0 and 1."
        );
    }

    if (albedoSurface < 0 || albedoSurface > 1) {
        throw new Error(
            "Surface albedo must be between 0 and 1."
        );
    }

    if (CO2 < 0 || CH4 < 0 || N2O < 0) {
        throw new Error(
            "Greenhouse-gas concentrations cannot be negative."
        );
    }

    if (albedoAtmosphere + SWabsAtmosphere > 1) {
        throw new Error(
            "Atmospheric albedo + atmospheric SW absorptivity cannot exceed 1."
        );
    }


    // =========================================================
    // REFERENCE SHORTWAVE BUDGET
    // =========================================================

    const SW_TOA_0 = S0 / 4;

    const SW_reflected_atmosphere_0 =
        aa0 * SW_TOA_0;

    const SW_absorbed_atmosphere_0 =
        SWabsAtmosphere * SW_TOA_0;

    const SW_surface_in_0 =
        (1 - aa0 - SWabsAtmosphere) * SW_TOA_0;

    const SW_reflected_surface_0 =
        as0 * SW_surface_in_0;

    const SW_absorbed_surface_0 =
        (1 - as0) * SW_surface_in_0;

    const ASR0 =
        SW_absorbed_atmosphere_0 +
        SW_absorbed_surface_0;


    // =========================================================
    // CURRENT SHORTWAVE BUDGET
    // =========================================================

    const SW_TOA =
        solarConstant / 4;

    const SW_reflected_atmosphere =
        albedoAtmosphere * SW_TOA;

    const SW_absorbed_atmosphere =
        SWabsAtmosphere * SW_TOA;

    const SW_surface_in =
        (1 - albedoAtmosphere - SWabsAtmosphere) *
        SW_TOA;

    const SW_reflected_surface =
        albedoSurface * SW_surface_in;

    const SW_absorbed_surface =
        (1 - albedoSurface) * SW_surface_in;

    const ASR =
        SW_absorbed_atmosphere +
        SW_absorbed_surface;


    // =========================================================
    // GREENHOUSE-GAS EFFECTIVE RADIATIVE FORCING
    // =========================================================

    const F_current =
        forcingGHGMeinshausen(CO2, CH4, N2O);

    const F_reference =
        forcingGHGMeinshausen(CO2_0, CH4_0, N2O_0);

    const dF_GHG =
        F_current - F_reference;


    // =========================================================
    // EQUILIBRIUM SURFACE TEMPERATURE
    // =========================================================

    const deltaT =
        (ASR - ASR0 + dF_GHG) / lambda;

    const Ts_K =
        T0_K + deltaT;

    const Ts_C =
        Ts_K - 273.15;


    // =========================================================
    // SURFACE LONGWAVE RADIATION
    // =========================================================

    const LW_surface =
        sigma * Math.pow(Ts_K, 4);


    // =========================================================
    // ATMOSPHERIC WINDOW
    // =========================================================

    const LW_surface_0 =
        sigma * Math.pow(T0_K, 4);

    const epsilonAbs =
        1 - LWwindow0 / LW_surface_0;

    const LW_window =
        (1 - epsilonAbs) * LW_surface;

    const LW_absorbed_atmosphere =
        epsilonAbs * LW_surface;


    // =========================================================
    // SENSIBLE AND LATENT HEAT
    // =========================================================

    let H =
        H0 * (1 + gammaH * deltaT);

    let LE =
        LE0 * (1 + gammaLE * deltaT);

    H = Math.max(H, 0);
    LE = Math.max(LE, 0);


    // =========================================================
    // OUTGOING LONGWAVE RADIATION
    // =========================================================

    const OLR =
        ASR0 + lambda * deltaT - dF_GHG;


    // =========================================================
    // ATMOSPHERIC DOWNWARD LW
    // =========================================================

    const LW_atmosphere_down =
        LW_surface +
        H +
        LE -
        SW_absorbed_surface;


    // =========================================================
    // ATMOSPHERIC UPWARD LW
    // =========================================================

    const LW_atmosphere_up =
        OLR - LW_window;


    // =========================================================
    // ENERGY-BALANCE CHECKS
    // =========================================================

    const TOA_imbalance =
        ASR - OLR;

    const surface_imbalance =
        SW_absorbed_surface +
        LW_atmosphere_down -
        LW_surface -
        H -
        LE;

    const atmosphere_imbalance =
        SW_absorbed_atmosphere +
        LW_absorbed_atmosphere +
        H +
        LE -
        LW_atmosphere_down -
        LW_atmosphere_up;


    // =========================================================
    // OUTPUT
    // =========================================================

    return {

        surfaceTemperature_C: Ts_C,
        surfaceTemperature_K: Ts_K,

        SW_TOA: SW_TOA,
        SW_reflected_atmosphere: SW_reflected_atmosphere,
        SW_absorbed_atmosphere: SW_absorbed_atmosphere,
        SW_reaching_surface: SW_surface_in,
        SW_reflected_surface: SW_reflected_surface,
        SW_absorbed_surface: SW_absorbed_surface,

        LW_surface: LW_surface,

        sensibleHeat: H,
        latentHeat: LE,

        LW_window: LW_window,
        LW_absorbed_atmosphere: LW_absorbed_atmosphere,

        LW_atmosphere_down: LW_atmosphere_down,
        LW_atmosphere_up: LW_atmosphere_up,

        ASR: ASR,
        OLR: OLR,

        GHGforcing_current: F_current,
        GHGforcing_reference: F_reference,
        dF_GHG: dF_GHG,

        TOA_imbalance: TOA_imbalance,
        surface_imbalance: surface_imbalance,
        atmosphere_imbalance: atmosphere_imbalance,

        epsilonAbs: epsilonAbs
    };
}

// =============================================================
// WEB INTERFACE
// =============================================================

const calculateButton =
    document.getElementById("calculateButton");

const resetButton =
    document.getElementById("resetButton");

const surfaceTemperatureDisplay =
    document.getElementById("surfaceTemperature");

// Flux displays
const fluxSWTOA =
    document.getElementById("fluxSWTOA");

const fluxSWReflectedAtmosphere =
    document.getElementById("fluxSWReflectedAtmosphere");

const fluxSWAbsorbedAtmosphere =
    document.getElementById("fluxSWAbsorbedAtmosphere");

const fluxSWReflectedSurface =
    document.getElementById("fluxSWReflectedSurface");

const fluxSWAbsorbedSurface =
    document.getElementById("fluxSWAbsorbedSurface");

const fluxLWSurface =
    document.getElementById("fluxLWSurface");

const fluxLWAbsorbedAtmosphere =
    document.getElementById("fluxLWAbsorbedAtmosphere");

const fluxLWAtmosphereDown =
    document.getElementById("fluxLWAtmosphereDown");

const fluxLWAtmosphereUp =
    document.getElementById("fluxLWAtmosphereUp");

const fluxLWWindow =
    document.getElementById("fluxLWWindow");

const fluxOLR =
    document.getElementById("fluxOLR");

const fluxSensibleHeat =
    document.getElementById("fluxSensibleHeat");

const fluxLatentHeat =
    document.getElementById("fluxLatentHeat");

// -------------------------------------------------------------
// Update all model outputs on the page
// -------------------------------------------------------------

function updateDisplay(result) {

    // Surface temperature
    surfaceTemperatureDisplay.textContent =
        result.surfaceTemperature_C.toFixed(1);


    // Shortwave fluxes
    fluxSWTOA.textContent =
        result.SW_TOA.toFixed(1);

    fluxSWReflectedAtmosphere.textContent =
        result.SW_reflected_atmosphere.toFixed(1);

    fluxSWAbsorbedAtmosphere.textContent =
        result.SW_absorbed_atmosphere.toFixed(1);

    fluxSWReflectedSurface.textContent =
        result.SW_reflected_surface.toFixed(1);

    fluxSWAbsorbedSurface.textContent =
        result.SW_absorbed_surface.toFixed(1);


    // Longwave fluxes
    fluxLWSurface.textContent =
        result.LW_surface.toFixed(1);

    fluxLWAbsorbedAtmosphere.textContent =
        result.LW_absorbed_atmosphere.toFixed(1);

    fluxLWAtmosphereDown.textContent =
        result.LW_atmosphere_down.toFixed(1);

    fluxLWAtmosphereUp.textContent =
        result.LW_atmosphere_up.toFixed(1);

    fluxLWWindow.textContent =
        result.LW_window.toFixed(1);

    fluxOLR.textContent =
        result.OLR.toFixed(1);


    // Turbulent fluxes
    fluxSensibleHeat.textContent =
        result.sensibleHeat.toFixed(1);

    fluxLatentHeat.textContent =
        result.latentHeat.toFixed(1);
}

// -------------------------------------------------------------
// Calculate button
// -------------------------------------------------------------

calculateButton.addEventListener("click", function () {

    // Read values from input fields
    const solarConstant =
        Number(document.getElementById("solarConstant").value);

    const albedoAtmosphere =
        Number(document.getElementById("albedoAtmosphere").value);

    const albedoSurface =
        Number(document.getElementById("albedoSurface").value);

    const CO2 =
        Number(document.getElementById("CO2").value);

    const CH4 =
        Number(document.getElementById("CH4").value);

    const N2O =
        Number(document.getElementById("N2O").value);

// Check web-app input ranges

if (solarConstant < 500 || solarConstant > 1500) {
    alert("Solar constant must be between 500 and 1500 W m⁻².");
    return;
}

if (albedoAtmosphere < 0 || albedoAtmosphere > 0.7) {
    alert("Atmospheric albedo must be between 0 and 0.7.");
    return;
}

if (albedoSurface < 0 || albedoSurface > 0.7) {
    alert("Surface albedo must be between 0 and 0.7.");
    return;
}

if (CO2 < 1 || CO2 > 30000) {
    alert("CO₂ concentration must be between 1 and 30,000 ppm.");
    return;
}

if (CH4 < 0 || CH4 > 10000) {
    alert("CH₄ concentration must be between 0 and 10,000 ppb.");
    return;
}

if (N2O < 0 || N2O > 2000) {
    alert("N₂O concentration must be between 0 and 2,000 ppb.");
    return;
}

    try {

        // Run radiation-balance model
        const result = radiationBalanceModel(
            solarConstant,
            albedoAtmosphere,
            albedoSurface,
            CO2,
            CH4,
            N2O
        );

        // Update all displayed model outputs
	updateDisplay(result);

    } catch (error) {

        alert(error.message);
    }

});

// -------------------------------------------------------------
// Reset button
// -------------------------------------------------------------

resetButton.addEventListener("click", function () {

    // Restore reference-climate input values
    document.getElementById("solarConstant").value = 1361;
    document.getElementById("albedoAtmosphere").value = "0.230";
    document.getElementById("albedoSurface").value = "0.120";
    document.getElementById("CO2").value = 424;
    document.getElementById("CH4").value = 1935;
    document.getElementById("N2O").value = 339;

    // Restore slider positions
    albedoAtmosphereSlider.value = 0.23;
    albedoSurfaceSlider.value = 0.120;


    // Calculate reference climate
    const result = radiationBalanceModel(
        1361,
        0.23,
        0.120,
        424,
        1935,
        339
    );

   // Update all displayed model outputs
      updateDisplay(result);

});

// -------------------------------------------------------------
// Atmospheric albedo slider
// -------------------------------------------------------------

const albedoAtmosphereSlider =
    document.getElementById("albedoAtmosphereSlider");

const albedoAtmosphereInput =
    document.getElementById("albedoAtmosphere");


// Slider -> numeric field
albedoAtmosphereSlider.addEventListener("input", function () {

    albedoAtmosphereInput.value =
        Number(albedoAtmosphereSlider.value).toFixed(3);

});


// Numeric field -> slider
albedoAtmosphereInput.addEventListener("input", function () {

    albedoAtmosphereSlider.value =
        albedoAtmosphereInput.value;

});


// Clamp manually entered value to allowed range
albedoAtmosphereInput.addEventListener("change", function () {

    let value = Number(albedoAtmosphereInput.value);

    // Clamp to allowed range
    value = Math.max(0, Math.min(0.7, value));

    // Update numeric field and slider
    albedoAtmosphereInput.value = value.toFixed(3);
    albedoAtmosphereSlider.value = value;

});

// -------------------------------------------------------------
// Surface albedo slider
// -------------------------------------------------------------

const albedoSurfaceSlider =
    document.getElementById("albedoSurfaceSlider");

const albedoSurfaceInput =
    document.getElementById("albedoSurface");


// Slider -> numeric field
albedoSurfaceSlider.addEventListener("input", function () {

    albedoSurfaceInput.value =
        Number(albedoSurfaceSlider.value).toFixed(3);

});


// Numeric field -> slider
albedoSurfaceInput.addEventListener("input", function () {

    albedoSurfaceSlider.value =
        albedoSurfaceInput.value;

});


// Clamp manually entered value to allowed range
albedoSurfaceInput.addEventListener("change", function () {

    let value = Number(albedoSurfaceInput.value);

    // Clamp to allowed range
    value = Math.max(0, Math.min(0.7, value));

    // Update numeric field and slider
    albedoSurfaceInput.value = value.toFixed(3);
    albedoSurfaceSlider.value = value;

});

// -------------------------------------------------------------
// Greenhouse-gas input fields: use integer values
// -------------------------------------------------------------

const CO2Input = document.getElementById("CO2");
const CH4Input = document.getElementById("CH4");
const N2OInput = document.getElementById("N2O");


CO2Input.addEventListener("change", function () {

    let value = Math.round(Number(CO2Input.value));

    value = Math.max(1, Math.min(30000, value));

    CO2Input.value = value;
});


CH4Input.addEventListener("change", function () {

    let value = Math.round(Number(CH4Input.value));

    value = Math.max(0, Math.min(10000, value));

    CH4Input.value = value;
});


N2OInput.addEventListener("change", function () {

    let value = Math.round(Number(N2OInput.value));

    value = Math.max(0, Math.min(2000, value));

    N2OInput.value = value;
});