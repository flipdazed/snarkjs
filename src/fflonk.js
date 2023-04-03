/*
    Copyright 2022 iden3 association.

    This file is part of snarkjs.

    snarkjs is a free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    snarkjs is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along with
    snarkjs. If not, see <https://www.gnu.org/licenses/>.
*/

// FFlonk constants
export const FF_T_POL_DEG_MIN = 3;

// FFlonk A
export const A = 12;

// ZKEY constants
export const ZKEY_FF_NSECTIONS = 17;

export const ZKEY_FF_HEADER_SECTION = 2;
export const ZKEY_FF_ADDITIONS_SECTION = 3;
export const ZKEY_FF_A_MAP_SECTION = 4;
export const ZKEY_FF_B_MAP_SECTION = 5;
export const ZKEY_FF_C_MAP_SECTION = 6;
export const ZKEY_FF_QL_SECTION = 7;
export const ZKEY_FF_QR_SECTION = 8;
export const ZKEY_FF_QM_SECTION = 9;
export const ZKEY_FF_QO_SECTION = 10;
export const ZKEY_FF_QC_SECTION = 11;
export const ZKEY_FF_SIGMA1_SECTION = 12;
export const ZKEY_FF_SIGMA2_SECTION = 13;
export const ZKEY_FF_SIGMA3_SECTION = 14;
export const ZKEY_FF_LAGRANGE_SECTION = 15;
export const ZKEY_FF_PTAU_SECTION = 16;
export const ZKEY_FF_C0_SECTION = 17;


export function computeLagrangeLiSi(roots, x, xi, curve) {
    const Fr = curve.Fr;
    const len = roots.length;

    const num = Fr.sub(Fr.exp(x, len), xi);
    const den1 = Fr.mul(Fr.e(len), Fr.exp(roots[0], len - 2));

    const Li = [];
    for (let i = 0; i < len; i++) {
        const den2 = roots[((len - 1) * i) % len];
        const den3 = Fr.sub(x, roots[i]);

        Li[i] = Fr.div(num, Fr.mul(Fr.mul(den1, den2), den3));
    }

    return Li;
}

export function computeLagrangeLiS2(roots, value, xi0, xi1, curve) {
    const Fr = curve.Fr;

    const Li = [];

    const len = roots[0].length;
    const n = len * roots.length;

    const num1 = Fr.exp(value, n);
    const num2 = Fr.mul(Fr.add(xi0, xi1), Fr.exp(value, len));
    const num3 = Fr.mul(xi0, xi1);
    const num = Fr.add(Fr.sub(num1, num2), num3);

    let den1 = Fr.mul(Fr.mul(Fr.e(len), roots[0][0]), Fr.sub(xi0, xi1));
    for (let i = 0; i < len; i++) {
        const den2 = roots[0][(len - 1) * i % len];
        const den3 = Fr.sub(value, roots[0][i]);

        const den = Fr.mul(den1,Fr.mul(den2, den3));

        Li[i] = Fr.div(num, den);
    }

    den1 = Fr.mul(Fr.mul(Fr.e(len), roots[1][0]), Fr.sub(xi1, xi0));
    for (let i = 0; i < len; i++) {
        const den2 = roots[1][(len - 1) * i % len];
        const den3 = Fr.sub(value, roots[1][i]);

        const den = Fr.mul(den1,Fr.mul(den2, den3));

        Li[i + len] = Fr.div(num, den);
    }

    return Li;
}
