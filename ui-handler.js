// UI-HANDLER.JS: FLOW CALCULATOR INPUTS, VALIDATION, AND OUTPUT DISPLAY

// Import all required functions from the flow-solver module
import {
    computeFlowProperties,
    solveMFromRatio,
    computeNormalShockRelations,
    solveM1FromShockRatio,
    computeObliqueShockProperties,
    solveBetaFromTurnAngle
} from './flow-solver.js';

const degToRad = Math.PI / 180; // Degree to radian conversion

// Clear output grid (DRY helper)
function clearGrid(gridSelector) {
    document.querySelectorAll(gridSelector + ' span').forEach(el => el.textContent = '');
}

document.addEventListener('DOMContentLoaded', () => {
    document.getElementById('calc-iso').addEventListener('click', calculateIsentropic);
    document.getElementById('calc-shock').addEventListener('click', calculateNormalShock);
    document.getElementById('calc-oblique').addEventListener('click', calculateObliqueShock);

    // Initial load calculations
    calculateIsentropic();
    calculateNormalShock();
    calculateObliqueShock();
});


// function to calculate Isentropic Flow Properties
function calculateIsentropic() {
    clearGrid('.results-grid:not(.shock-grid):not(.oblique-grid)');

    const inputType = document.getElementById('input-select-iso').value;
    const inputValue = parseFloat(document.getElementById('input-value-iso').value);
    const gamma = parseFloat(document.getElementById('gamma-input-iso').value);

    // Input validation: Alert on NaN
    if (isNaN(inputValue) || isNaN(gamma)) {
        window.alert("Error: Invalid or missing Input value or Gamma.");
        document.getElementById('M_output_iso').textContent = 'Error';
        return;
    }

    /* * START OF NEW ISENTROPIC VALIDATION 
    */
    let error_message = null;
    if (gamma <= 1.0) {
        error_message = "Error: Gamma (ratio of specific heats) must be greater than 1.0.";
    } else if (inputType.endsWith('_sup') && inputValue < 1.0) {
         // Check if A/A* is chosen as supersonic but value is < 1
         error_message = "Error: Area ratio (A/A*) for the supersonic branch must be $\\geq 1.0$ (Critical area).";
    } else if (inputType === 'T_T0' && (inputValue < 0 || inputValue > 1)) {
        error_message = "Error: T/T₀ must be between 0 and 1 (Isentropic relations).";
    } else if (inputType === 'p_p0' && (inputValue < 0 || inputValue > 1)) {
        error_message = "Error: p/p₀ must be between 0 and 1 (Isentropic relations).";
    } else if (inputType === 'rho_rho0' && (inputValue < 0 || inputValue > 1)) {
        error_message = "Error: $\\rho/\\rho₀$ must be between 0 and 1 (Isentropic relations).";
    } else if (inputType === 'mu_deg' && (inputValue <= 0 || inputValue > 90)) {
        error_message = "Error: Mach angle ($\\mu$) must be between 0° and 90° for Mach numbers $\\geq 1$.";
    } else if (inputType === 'nu_deg' && inputValue < 0) {
        // Upper limit is very large (e.g., ~130 for gamma=1.4)
        error_message = "Error: Prandtl-Meyer function ($\\nu$) must be $\\geq 0$.";
    }

    if (error_message) {
        window.alert(error_message);
        document.getElementById('M_output_iso').textContent = 'Error';
        return;
    }
    /*
    * END OF NEW ISENTROPIC VALIDATION
    */

    let M;
    if (inputType === 'M') {
        M = inputValue;
    } else {
        M = solveMFromRatio(inputType, inputValue, gamma);
    }
    
    // Check for invalid M solution: Alert on NaN or M < 0
    if (isNaN(M) || M < 0) {
        /* * IMPROVED ERROR MESSAGE FOR UN-PHYSICAL RATIO VALUE
        */
        window.alert(`Error: No valid Mach number (M) found for the given ratio (${inputType} = ${inputValue}). Check the physical range for the ratio.`);
        document.getElementById('M_output_iso').textContent = 'Error';
        return;
    }

    // Compute all isentropic flow properties
    const results = computeFlowProperties(M, gamma);

    // Display Results for Isentropic
    document.getElementById('M_output_iso').textContent = results.M.toFixed(6);
    // Handle mu (Mach angle) display for subsonic flow
    document.getElementById('mu_deg_output_iso').textContent = (results.mu_deg !== 'N/A (<M1)') ? results.mu_deg.toFixed(6) : 'N/A (<M₁)';
    document.getElementById('nu_deg_output_iso').textContent = results.nu_deg.toFixed(6);
    document.getElementById('T_T0_output_iso').textContent = results.T_T0.toFixed(6);
    document.getElementById('p_p0_output_iso').textContent = results.p_p0.toFixed(6);
    document.getElementById('rho_rho0_output_iso').textContent = results.rho_rho0.toFixed(6);
    document.getElementById('T_Tstar_output_iso').textContent = results.T_Tstar.toFixed(6);
    document.getElementById('p_pstar_output_iso').textContent = results.p_pstar.toFixed(6);
    document.getElementById('rho_rhostar_output_iso').textContent = results.rho_rhostar.toFixed(6);
    document.getElementById('A_Astar_output_iso').textContent = results.A_Astar.toFixed(6);
}


// Function to calculate Normal Shock Relations
function calculateNormalShock() {
    // Clear normal shock grid
    clearGrid('.shock-grid');

    const inputType = document.getElementById('input-select-shock').value;
    const inputValue = parseFloat(document.getElementById('input-value-shock').value);
    const gamma = parseFloat(document.getElementById('gamma-input-shock').value);

    // Input validation: Alert on NaN
    if (isNaN(inputValue) || isNaN(gamma)) {
        window.alert("Error: Invalid or missing Input value or Gamma.");
        document.getElementById('M1_output_shock').textContent = 'Error';
        return;
    }

    /* * START OF NEW NORMAL SHOCK VALIDATION 
    */
    let shock_error = null;
    if (gamma <= 1.0) {
        shock_error = "Error: Gamma must be greater than 1.0.";
    } else if (inputType === 'p2_p1' && inputValue < 1) {
        shock_error = "Error: The static pressure ratio ($p₂/p₁$) across a shock must be $\\geq 1.0$.";
    } else if (inputType === 'rho2_rho1' && (inputValue < 1 || inputValue > (gamma + 1) / (gamma - 1) + 1e-6)) {
        shock_error = `Error: The density ratio ($\\rho₂/\\rho₁$) must be $\\geq 1.0$ and less than the maximum possible limit of $${((gamma + 1) / (gamma - 1)).toFixed(3)}$.`;
    } else if (inputType === 'T2_T1' && inputValue < 1) {
        shock_error = "Error: The temperature ratio ($T₂/T₁$) must be $\\geq 1.0$.";
    } else if (inputType === 'p02_p01' && (inputValue > 1 || inputValue <= 0)) {
        shock_error = "Error: The stagnation pressure ratio ($p₀₂/p₀₁$) must be $0 < p₀₂/p₀₁ \\leq 1$ (Losses).";
    } else if (inputType === 'p1_p02' && (inputValue < 0 || inputValue > 1)) {
        shock_error = "Error: $p₁/p₀₂$ must be between 0 and 1.";
    }

    if (shock_error) {
        window.alert(shock_error);
        document.getElementById('M1_output_shock').textContent = 'Error';
        return;
    }
    /* * END OF NEW NORMAL SHOCK VALIDATION 
    */

    let M1;

    // 1. Determine M1 based on input type
    M1 = (inputType === 'M1') ? inputValue : solveM1FromShockRatio(inputType, inputValue, gamma);
    
    // Check if the inversion failed (e.g., for M2, M1 is complex)
    if (isNaN(M1)) {
        window.alert(`Error: No valid M₁ could be found for the given ratio (${inputType} = ${inputValue}). The value may be outside the physical range for a normal shock.`);
        document.getElementById('M1_output_shock').textContent = 'Error';
        return;
    }

    // 2. Handle Subsonic/Invalid M1 (no normal shock)
    if (M1 < 1.0001) {
        
        const isM1Direct = (inputType === 'M1');
        const M1_display = isM1Direct ? inputValue : M1;

        if (isM1Direct || inputType !== 'M2') { // Alert if M1 is directly input as subsonic, or if ratio inversion landed on M1~1
            // Use the input value for the computation to get the correct p/p0
            const results = computeNormalShockRelations(inputValue, gamma);

            // Display Subsonic/Isentropic Results (ratios=1)
            document.getElementById('M1_output_shock').textContent = M1_display.toFixed(6);
            document.getElementById('M2_output_shock').textContent = 'N/A (No Shock)'; // Clearer UI
            document.getElementById('p2_p1_output_shock').textContent = '1.0000';
            document.getElementById('rho2_rho1_output_shock').textContent = '1.0000';
            document.getElementById('T2_T1_output_shock').textContent = '1.0000';
            document.getElementById('p02_p01_output_shock').textContent = '1.0000';
            document.getElementById('p1_p02_output_shock').textContent = results.p1_p02.toFixed(6);
            
            // Alert for subsonic or invalid M1
            window.alert(`M₁=${M1_display.toFixed(6)} is subsonic. A normal shock cannot form. Results reflect isentropic flow (ratios = 1).`);
            return;
        }
    }

    // 3. Compute All Properties (Supersonic M1)
    const results = computeNormalShockRelations(M1, gamma);

    // 4. Display Supersonic Results
    document.getElementById('M1_output_shock').textContent = results.M1.toFixed(6);
    document.getElementById('M2_output_shock').textContent = results.M2.toFixed(6);
    document.getElementById('p2_p1_output_shock').textContent = results.p2_p1.toFixed(6);
    document.getElementById('rho2_rho1_output_shock').textContent = results.rho2_rho1.toFixed(6);
    document.getElementById('T2_T1_output_shock').textContent = results.T2_T1.toFixed(6);
    document.getElementById('p02_p01_output_shock').textContent = results.p02_p01.toFixed(6);
    document.getElementById('p1_p02_output_shock').textContent = results.p1_p02.toFixed(6);
}


// Function to calculate Oblique Shock Relations (NO CHANGE NEEDED HERE)
function calculateObliqueShock() {
    clearGrid('.oblique-grid');

    const inputType = document.getElementById('input-select-oblique').value;
    const M1 = parseFloat(document.getElementById('input-M1-oblique').value);
    const inputValue = parseFloat(document.getElementById('input-value-oblique').value);
    const gamma = parseFloat(document.getElementById('gamma-input-oblique').value);

    // Input validation: Alert on NaN or subsonic M1
    if (isNaN(M1) || M1 < 1.0001 || isNaN(inputValue) || isNaN(gamma)) {
        window.alert("Error: Invalid Mach number (M₁) or missing Input value/Gamma.");
        document.getElementById('M2_output_oblique').textContent = 'Error';
        return;
    }
    
    let beta_deg;
    let turn_angle_deg;
    let results;
    let M1n_calc = NaN;

    if (inputType === 'wave_angle') {
        beta_deg = inputValue;
        // Check if beta is in the valid range for an oblique shock
        const mu_deg = Math.asin(1 / M1) * (180 / Math.PI);
        if (beta_deg <= mu_deg || beta_deg > 90) {
            window.alert(`Error: Wave angle ($\\beta$) must be between the Mach angle ($\\,${mu_deg.toFixed(2)}^\\circ$) and $90^\\circ$.`);
            document.getElementById('M2_output_oblique').textContent = 'Error';
            return;
        }

        // Compute properties from M1 and beta
        results = computeObliqueShockProperties(M1, beta_deg, gamma);
        M1n_calc = M1 * Math.sin(beta_deg * degToRad);
        
    } else if (inputType.startsWith('turn_angle')) {
        turn_angle_deg = inputValue;
        
        // Solve for wave angles (beta) from M1 and turn angle (theta)
        const [beta_weak, beta_strong] = solveBetaFromTurnAngle(M1, turn_angle_deg, gamma);
        
        if (isNaN(beta_weak) && isNaN(beta_strong)) {
            // Detached shock error: Alert user
            window.alert("Error: Detached Shock. The flow deflection angle ($\\theta$) is too high for M₁.");
            document.getElementById('M2_output_oblique').textContent = 'Detached';
            document.getElementById('turn_angle_output_oblique').textContent = turn_angle_deg.toFixed(6);
            return;
        }

        // Select weak or strong solution
        beta_deg = (inputType === 'turn_angle_weak') ? beta_weak : beta_strong;
        
        if (isNaN(beta_deg)) {
            // Selected solution doesn't exist (e.g., attempting to find weak solution at theta_max)
            document.getElementById('M2_output_oblique').textContent = (inputType === 'turn_angle_weak') ?
                'N/A (No Weak)' : 'N/A (No Strong)';
            document.getElementById('turn_angle_output_oblique').textContent = turn_angle_deg.toFixed(6);
            return;
        }
        
        // Compute properties for the chosen beta
        results = computeObliqueShockProperties(M1, beta_deg, gamma);
        M1n_calc = M1 * Math.sin(beta_deg * degToRad);
        
    } else if (inputType === 'M1n') {
        M1n_calc = inputValue;
        // M1n must be >= 1 and <= M1: Alert on invalid input
        if (M1n_calc > M1) {
            window.alert("Error: M₁n cannot be greater than M₁ (M₁n must be $\\leq$ M₁).");
            document.getElementById('M2_output_oblique').textContent = 'Error';
            return;
        }
        if (M1n_calc < 1.0001) {
            window.alert("Error: M₁n must be supersonic ($\\,> 1.0$) for an oblique shock.");
            document.getElementById('M2_output_oblique').textContent = 'Error';
            return;
        }
        // Calculate beta from M1n and M1
        beta_deg = Math.asin(M1n_calc / M1) * (180 / Math.PI);
        results = computeObliqueShockProperties(M1, beta_deg, gamma);
    }

    // Always display M1n first
    if (!isNaN(M1n_calc)) {
        document.getElementById('M1n_output_oblique').textContent = M1n_calc.toFixed(4);
    }

    // Check if the oblique shock calculation failed internally (should be caught by checks above, but good safeguard)
    if (!results || !results.valid) {
        window.alert(`Error: Invalid Oblique Shock. M₁n must be supersonic ($\\,> 1.0$) or another detachment condition has occurred. Details: ${results.error || 'Unknown'}`);
        document.getElementById('M2_output_oblique').textContent = 'Error';
        return;
    }
    
    // Display Results
    document.getElementById('M2_output_oblique').textContent = results.M2.toFixed(6);
    document.getElementById('turn_angle_output_oblique').textContent = results.turn_angle_deg.toFixed(6);
    document.getElementById('wave_angle_output_oblique').textContent = results.wave_angle_deg.toFixed(6);
    document.getElementById('M2n_output_oblique').textContent = results.M2n.toFixed(6);

    document.getElementById('p2_p1_output_oblique').textContent = results.p2_p1.toFixed(6);
    document.getElementById('rho2_rho1_output_oblique').textContent = results.rho2_rho1.toFixed(6);
    document.getElementById('T2_T1_output_oblique').textContent = results.T2_T1.toFixed(6);
    document.getElementById('p02_p01_output_oblique').textContent = results.p02_p01.toFixed(6);
}