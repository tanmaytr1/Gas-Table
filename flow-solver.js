// FLOW-SOLVER.JS: CORE COMPRESSIBLE FLOW MATHEMATICS

const maxIter = 100; // Max iterations for numerical solutions
const tolerance = 1e-6; // Tolerance for numerical solutions

const degToRad = Math.PI / 180; // Degrees to radians
const radToDeg = 180 / Math.PI; // Radians to degrees

// ISENTROPIC FLOW RELATIONS

/**
 * Calculates all isentropic flow properties (ratios to stagnation and critical states) for M and gamma.
 */
export function computeFlowProperties(M, gamma) {
    const R_g = (gamma - 1) / 2; 
    const term1 = 1 + R_g * M * M; 

    // Critical relations (M=1)
    const T_star_T0 = 2 / (gamma + 1);
    const p_star_p0 = Math.pow(T_star_T0, gamma / (gamma - 1));
    const rho_star_rho0 = Math.pow(T_star_T0, 1 / (gamma - 1));

    // Stagnation relations (M=M)
    const T_T0 = 1 / term1;
    const p_p0 = Math.pow(T_T0, gamma / (gamma - 1));
    const rho_rho0 = Math.pow(T_T0, 1 / (gamma - 1));

    // Ratios to critical state
    const T_Tstar = T_T0 / T_star_T0;
    const p_pstar = p_p0 / p_star_p0;
    const rho_rhostar = rho_rho0 / rho_star_rho0;

    // Area ratio (A/A*)
    const A_Astar = (1 / M) * Math.pow(
        T_star_T0 * term1,
        (gamma + 1) / (2 * (gamma - 1))
    );

    let mu_deg = 'N/A (<M1)'; // Mach angle (mu)
    let nu_deg = 0.000; // Prandtl-Meyer function (nu)

    if (M >= 1.0) {
        mu_deg = Math.asin(1 / M) * radToDeg; // mu = arcsin(1/M)
        const G_r = (gamma + 1) / (gamma - 1);
        const pm_term_sq = M * M - 1;
        const K_pm = Math.sqrt(G_r);
        const arg1 = Math.sqrt((gamma - 1) / (gamma + 1) * pm_term_sq);
        const arg2 = Math.sqrt(pm_term_sq);
        // Prandtl-Meyer function calculation in radians
        const nu_rad = K_pm * Math.atan(arg1) - Math.atan(arg2);
        nu_deg = nu_rad * radToDeg;
    }

    return { M, T_T0, p_p0, rho_rho0, T_Tstar, p_pstar, rho_rhostar, A_Astar, mu_deg, nu_deg };
}

/**
 * Solves for Mach number M given a specific isentropic ratio.
 */
export function solveMFromRatio(ratio_name, ratio_value, gamma) {
    const R_g = (gamma - 1) / 2;
    let M_squared;

    // Numerical solutions for non-linear relations
    switch (ratio_name) {
        case 'T_T0': M_squared = (1 / ratio_value - 1) / R_g; break;
        case 'p_p0': const E_p = gamma / (gamma - 1); M_squared = (Math.pow(ratio_value, -1 / E_p) - 1) / R_g; break;
        case 'rho_rho0': const E_rho = 1 / (gamma - 1); M_squared = (Math.pow(ratio_value, -1 / E_rho) - 1) / R_g; break;
        case 'mu_deg': const mu_rad = ratio_value * degToRad; return 1 / Math.sin(mu_rad); // M = 1/sin(mu)
         
        case 'nu_deg': // Prandtl-Meyer inversion (Bisection method)
            const nu_target_rad = ratio_value * degToRad;
            let M_low = 1.0, M_high = 100.0;
            const pm_function = (M_in) => {
                const G_r = (gamma + 1) / (gamma - 1);
                const pm_term_sq = M_in * M_in - 1;
                const K_pm = Math.sqrt(G_r);
                const arg1 = Math.sqrt((gamma - 1) / (gamma + 1) * pm_term_sq);
                const arg2 = Math.sqrt(pm_term_sq);
                return K_pm * Math.atan(arg1) - Math.atan(arg2);
            };
            for(let i = 0; i < maxIter; i++) {
                const M_mid = (M_low + M_high) / 2;
                const nu_mid = pm_function(M_mid);
                if (Math.abs(nu_mid - nu_target_rad) < tolerance) return M_mid;
                if (nu_mid < nu_target_rad) M_low = M_mid;
                else M_high = M_mid;
            }
            return (M_low + M_high) / 2;
         
        case 'A_Astar_sub': 
        case 'A_Astar_sup': 
            // A/A* inversion (Simple search/Refined search, as it's non-monotonic)
            const isSupersonic = ratio_name.endsWith('_sup');
            let M_start = isSupersonic ? 1.01 : 0.01;
            let M_end = isSupersonic ? 10.0 : 0.99;
             
            const area_function = (M_in) => {
                const term1_in = 1 + R_g * M_in * M_in;
                const T_star_T0 = 2 / (gamma + 1);
                return (1 / M_in) * Math.pow(T_star_T0 * term1_in, (gamma + 1) / (2 * (gamma - 1)));
            };
            let M_guess = M_start;
            const step = (M_end - M_start) / maxIter; // Simple linear step search
            let min_diff = Infinity;
            let best_M = NaN;

            for (let i = 0; i < maxIter; i++) {
                const A_Astar_calc = area_function(M_guess);
                const diff = Math.abs(A_Astar_calc - ratio_value);
                if (diff < tolerance) return M_guess;
                if (diff < min_diff) {
                    min_diff = diff;
                    best_M = M_guess;
                }
                M_guess += step;
            }
            return best_M;
    }
     
    return (M_squared >= 0) ? Math.sqrt(M_squared) : NaN;
}



// NORMAL SHOCK RELATIONS


/**
 * Calculates all normal shock relations for M1 and gamma.
 */
export function computeNormalShockRelations(M1, gamma) {
    const M1_sq = M1 * M1;
    const gamma_minus_1 = gamma - 1;
    const gamma_plus_1 = gamma + 1;
    const two_gamma = 2 * gamma;
    const R_g = gamma_minus_1 / 2;

    // Subsonic M1 -> isentropic results
    if (M1 < 1.0001) {
          const iso = computeFlowProperties(M1, gamma);
          return {
              M1: M1, M2: M1, p2_p1: 1.0, rho2_rho1: 1.0, T2_T1: 1.0, p02_p01: 1.0, p1_p02: iso.p_p0
          };
    }

    // M₂ (Downstream Mach number)
    const M2_sq = (M1_sq + 2 / gamma_minus_1) / ((two_gamma / gamma_minus_1) * M1_sq - 1);
    const M2 = Math.sqrt(M2_sq);

    // p₂/p₁ (Static Pressure Ratio)
    const p2_p1 = 1 + (two_gamma / gamma_plus_1) * (M1_sq - 1);

    // ρ₂/ρ₁ (Density Ratio)
    const rho2_rho1 = (gamma_plus_1 * M1_sq) / (2 + gamma_minus_1 * M1_sq);

    // T₂/T₁ (Temperature Ratio)
    const T2_T1 = p2_p1 / rho2_rho1;

    // p₀₂/p₀₁ (Stagnation Pressure Ratio)
    const term_A = (gamma_plus_1 * M1_sq) / (2 + gamma_minus_1 * M1_sq);
    const term_B = (gamma_plus_1) / (two_gamma * M1_sq - gamma_minus_1);
    const p02_p01 = Math.pow(term_A, gamma / gamma_minus_1) * Math.pow(term_B, 1 / gamma_minus_1);

    // p₁/p₀₂ (Isentropic p₁/p₀₁ divided by p₀₂/p₀₁)
    const p1_p01 = Math.pow(1 + R_g * M1_sq, -gamma / gamma_minus_1);
    const p1_p02 = p1_p01 / p02_p01;

    return { M1, M2, p2_p1, rho2_rho1, T2_T1, p02_p01, p1_p02 };
}


/**
 * Solves for upstream Mach number M1 given a shock ratio.
 */

export function solveM1FromShockRatio(ratio_name, ratio_value, gamma) {
    const M_sq = ratio_value * ratio_value;
    const gamma_minus_1 = gamma - 1;
    const gamma_plus_1 = gamma + 1;
    const two_gamma = 2 * gamma;
    const R_g = gamma_minus_1 / 2;

    switch (ratio_name) {
        case 'M2':
            // M1^2 from M2^2
            const M1_sq_M2_num = 1 + R_g * M_sq;
            const M1_sq_M2_den = gamma * M_sq - R_g;
            const M1_sq_M2 = M1_sq_M2_num / M1_sq_M2_den;
            return (M1_sq_M2 >= 1) ? Math.sqrt(M1_sq_M2) : NaN;

        case 'p2_p1':
            // M1^2 from p2/p1
            const M1_sq_p = (ratio_value - 1) * (gamma_plus_1 / two_gamma) + 1;
            return (M1_sq_p >= 1) ? Math.sqrt(M1_sq_p) : NaN;

        case 'rho2_rho1':
            // M1^2 from rho2/rho1
            const M1_sq_rho_num = two_gamma * ratio_value;
            const M1_sq_rho_den = gamma_plus_1 - gamma_minus_1 * ratio_value;
            const M1_sq_rho = M1_sq_rho_num / M1_sq_rho_den;
            return (M1_sq_rho >= 1) ? Math.sqrt(M1_sq_rho) : NaN;

        case 'T2_T1':
            // T2/T1 inversion (Bisection method)
            let M1_low = 1.0;
            let M1_high = 100.0;

            const T_ratio_func = (M1_in) => {
                const M1_in_sq = M1_in * M1_in;
                const p2_p1_calc = 1 + (two_gamma / gamma_plus_1) * (M1_in_sq - 1);
                const rho2_rho1_calc = (gamma_plus_1 * M1_in_sq) / (2 + gamma_minus_1 * M1_in_sq);
                return p2_p1_calc / rho2_rho1_calc;
            };

            for(let i = 0; i < maxIter; i++) {
                const M1_mid = (M1_low + M1_high) / 2;
                const T_mid = T_ratio_func(M1_mid);
                if (Math.abs(T_mid - ratio_value) < tolerance) return M1_mid;
                if (T_mid < ratio_value) M1_low = M1_mid; // T ratio increases with M1
                else M1_high = M1_mid;
            }
            return (M1_low + M1_high) / 2;

        case 'p02_p01':
        case 'p1_p02':
            // Stagnation pressure ratio inversion (Bisection method)
            let M1_low_p0 = 1.0;
            let M1_high_p0 = 500.0;

            const P0_ratio_func = (M1_in) => {
                const results = computeNormalShockRelations(M1_in, gamma);
                return (ratio_name === 'p02_p01') ? results.p02_p01 : results.p1_p02;
            };

            for(let i = 0; i < maxIter; i++) {
                const M1_mid = (M1_low_p0 + M1_high_p0) / 2;
                const P0_mid = P0_ratio_func(M1_mid);

                let target = ratio_value;

                if (Math.abs(P0_mid - target) < tolerance) return M1_mid;

                // p02/p01 and p1/p02 decrease with M1, so logic is reversed
                if (P0_mid > target) M1_low_p0 = M1_mid;
                else M1_high_p0 = M1_mid;
            }
            return (M1_low_p0 + M1_high_p0) / 2;
    }

    return NaN;
}



// OBLIQUE SHOCK RELATIONS

// Calculates turn angle (theta) from M1 and wave angle (beta). (theta-beta-M relation)
function thetaBetaM(M1, beta_rad, gamma) {
    const M1_sq = M1 * M1;
    const sin_beta_sq = Math.sin(beta_rad) * Math.sin(beta_rad);

    // tan(theta) formula
    const tan_theta_num = 2 * (M1_sq * sin_beta_sq - 1) * (1 / Math.tan(beta_rad));
    const tan_theta_den = M1_sq * (gamma + Math.cos(2 * beta_rad)) + 2;

    return Math.atan(tan_theta_num / tan_theta_den);
}


// Calculates flow properties behind an oblique shock using M1 and wave angle (beta).
export function computeObliqueShockProperties(M1, beta_deg, gamma) {
    const beta_rad = beta_deg * degToRad;
    const sin_beta = Math.sin(beta_rad);
    const M1_sq = M1 * M1;

    // 1. Normal Mach Number Components
    const M1n = M1 * sin_beta;
    const M1n_sq = M1n * M1n;

    if (M1n < 1.0001) {
        // Shock is degenerate (Mach wave) or impossible
        return {
            valid: false, error: "M₁n is subsonic."
        };
    }

    // Constants for normal shock relations
    const gamma_minus_1 = gamma - 1;
    const two_gamma = 2 * gamma;
    const gamma_plus_1 = gamma + 1;

    // 2. Ratios based on M1n (Normal Shock Relations)
    const p2_p1 = 1 + (two_gamma / gamma_plus_1) * (M1n_sq - 1);
    const rho2_rho1 = (gamma_plus_1 * M1n_sq) / (2 + gamma_minus_1 * M1n_sq);
    const T2_T1 = p2_p1 / rho2_rho1;
    const term_A = (gamma_plus_1 * M1n_sq) / (2 + gamma_minus_1 * M1n_sq);
    const term_B = (gamma_plus_1) / (two_gamma * M1n_sq - gamma_minus_1);
    const p02_p01 = Math.pow(term_A, gamma / gamma_minus_1) * Math.pow(term_B, 1 / gamma_minus_1);

    // 3. Downstream Normal Mach Number
    const M2n_sq = (M1n_sq + 2 / gamma_minus_1) / ((two_gamma / gamma_minus_1) * M1n_sq - 1);
    const M2n = Math.sqrt(M2n_sq);

    // 4. Turn Angle (theta)
    const theta_rad = thetaBetaM(M1, beta_rad, gamma);
    const turn_angle_deg = theta_rad * radToDeg;

    if (theta_rad < 0) { // Indicates an error or mathematically impossible solution (e.g. M1n too high)
        return { valid: false, error: "Detached or Impossible Shock" };
    }

    // 5. Downstream Mach Number (M₂ = M₂n / sin(β - θ))
    const beta_minus_theta = beta_rad - theta_rad;

    if (Math.sin(beta_minus_theta) === 0) {
        return { valid: false, error: "Sine of beta-theta is zero." };
    }

    const M2 = M2n / Math.sin(beta_minus_theta);

    return {
        valid: true,
        M1: M1, M2: M2, turn_angle_deg, wave_angle_deg: beta_deg,
        M1n: M1n, M2n: M2n, p2_p1, rho2_rho1, T2_T1, p02_p01
    };
}


// Solves for the wave angle (beta) given M1 and the turn angle (theta). (Inversion of theta-beta-M)
export function solveBetaFromTurnAngle(M1, theta_deg, gamma) {
    if (M1 < 1.0001) return [NaN, NaN]; // No shock for subsonic flow

    const theta_rad = theta_deg * degToRad;
    const mu_rad = Math.asin(1 / M1); // Mach angle (minimum beta)
    const mu_deg = mu_rad * radToDeg;

    // Find the maximum turn angle (theta_max) and corresponding beta
    let beta_max_low = mu_deg, beta_max_high = 90.0;
    let beta_at_theta_max = mu_deg;
    let theta_max = 0.0;
     
    // Numerical search for theta_max (simple optimization search)
    for (let i = 0; i < maxIter; i++) {
        const beta_mid = (beta_max_low + beta_max_high) / 2;
        const theta_mid = thetaBetaM(M1, beta_mid * degToRad, gamma);
         
        // If theta_mid is greater than theta_max, the maximum must be at a higher beta
        if (theta_mid > theta_max) {
            theta_max = theta_mid;
            beta_at_theta_max = beta_mid;
            beta_max_low = beta_mid;
        } else {
            // The maximum has been passed
            beta_max_high = beta_mid;
        }
    }
    theta_max *= radToDeg;

    if (theta_deg > theta_max + tolerance) { // Check for detached shock
        return [NaN, NaN]; // Detached shock
    }

    let beta_weak = NaN;
    let beta_strong = NaN;

    // 1. Weak Shock Solution: Beta is between mu and beta_at_theta_max (Bisection)
    let low = mu_deg;
    let high = beta_at_theta_max;

    for (let i = 0; i < maxIter; i++) {
        const mid = (low + high) / 2;
        const theta_mid = thetaBetaM(M1, mid * degToRad, gamma) * radToDeg;
         
        if (Math.abs(theta_mid - theta_deg) < tolerance) {
            beta_weak = mid;
            break;
        }
        if (theta_mid > theta_deg) high = mid; // theta decreases as beta increases in this range
        else low = mid;
    }

    // 2. Strong Shock Solution: Beta is between beta_at_theta_max and 90 degrees (Bisection)
    low = beta_at_theta_max;
    high = 90.0;
     
    for (let i = 0; i < maxIter; i++) {
        const mid = (low + high) / 2;
        const theta_mid = thetaBetaM(M1, mid * degToRad, gamma) * radToDeg;
         
        if (Math.abs(theta_mid - theta_deg) < tolerance) {
            beta_strong = mid;
            break;
        }
        if (theta_mid < theta_deg) high = mid; // theta decreases as beta increases in this range
        else low = mid;
    }

    return [beta_weak, beta_strong];
}