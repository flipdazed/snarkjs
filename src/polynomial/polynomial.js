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

import {BigBuffer} from "ffjavascript";

export class Polynomial {
    constructor(coefficients = new Uint8Array(0), Fr, logger) {
        this.coef = coefficients;
        this.Fr = Fr;
        this.logger = logger;
    }

    static async fromBuffer(buffer, Fr, logger) {
        let coefficients = await Fr.ifft(buffer);

        return new Polynomial(coefficients, Fr, logger);
    }

    blindCoefficients(blindingFactors) {
        blindingFactors = blindingFactors || [];

        const blindedCoefficients = new BigBuffer((this.length() + blindingFactors.length) * this.Fr.n8);
        blindedCoefficients.set(this.coef, 0);
        for (let i = 0; i < blindingFactors.length; i++) {
            blindedCoefficients.set(
                this.Fr.add(
                    blindedCoefficients.slice((this.length() + i) * this.Fr.n8, (this.length() + i + 1) * this.Fr.n8),
                    blindingFactors[i]
                ),
                (this.length() + i) * this.Fr.n8
            );
            blindedCoefficients.set(
                this.Fr.sub(
                    blindedCoefficients.slice(i * this.Fr.n8, (i + 1) * this.Fr.n8),
                    blindingFactors[i]
                ),
                i * this.Fr.n8
            );
        }
        this.coef = blindedCoefficients;
    }

    getCoef(index) {
        const i_n8 = index * this.Fr.n8;

        if (i_n8 + this.Fr.n8 > this.coef.length) return this.Fr.zero;

        return this.coef.slice(i_n8, i_n8 + this.Fr.n8);
    }

    setCoef(index, value) {
        if (index > this.degree()) {
            throw new Error("Coef index is not available");
        }

        this.coef.set(value, index * this.Fr.n8);
    }

    static async to4T(buffer, domainSize, blindingFactors, Fr) {
        blindingFactors = blindingFactors || [];
        let a = await Fr.ifft(buffer);

        const a4 = new BigBuffer(domainSize * 4 * Fr.n8);
        a4.set(a, 0);

        const A4 = await Fr.fft(a4);

        if (blindingFactors.length === 0) {
            return [a, A4];
        }

        const a1 = new BigBuffer((domainSize + blindingFactors.length) * Fr.n8);
        a1.set(a, 0);
        for (let i = 0; i < blindingFactors.length; i++) {
            a1.set(
                Fr.add(
                    a1.slice((domainSize + i) * Fr.n8, (domainSize + i + 1) * Fr.n8),
                    blindingFactors[i]
                ),
                (domainSize + i) * Fr.n8
            );
            a1.set(
                Fr.sub(
                    a1.slice(i * Fr.n8, (i + 1) * Fr.n8),
                    blindingFactors[i]
                ),
                i * Fr.n8
            );
        }

        return [a1, A4];
    }

    length() {
        let length = this.coef.byteLength / this.Fr.n8;
        if (length !== Math.floor(this.coef.byteLength / this.Fr.n8)) {
            throw new Error("Polynomial coefficients buffer has incorrect size");
        }
        if (0 === length) {
            if (this.logger) {
                this.logger.warn("Polynomial has length zero");
            }
        }
        return length;
    }

    degree() {
        for (let i = this.length() - 1; i > 0; i--) {
            const i_n8 = i * this.Fr.n8;
            if (!this.Fr.eq(this.Fr.zero, this.coef.slice(i_n8, i_n8 + this.Fr.n8))) {
                return i;
            }
        }

        return 0;
    }

    evaluate(point) {
        let res = this.Fr.zero;

        for (let i = this.length(); i > 0; i--) {
            let i_n8 = (i - 1) * this.Fr.n8;
            const currentCoefficient = this.coef.slice(i_n8, i_n8 + this.Fr.n8);
            res = this.Fr.add(currentCoefficient, this.Fr.mul(res, point));
        }

        return res;
    }

    add(polynomial, blindingValue) {
        let other = false;

        if (polynomial.length() > this.length()) {
            other = true;
        }

        const thisLength = this.length();
        const polyLength = polynomial.length();
        for (let i = 0; i < Math.max(thisLength, polyLength); i++) {
            const i_n8 = i * this.Fr.n8;

            const a = i < thisLength ? this.coef.slice(i_n8, i_n8 + this.Fr.n8) : this.Fr.zero;
            let b = i < polyLength ? polynomial.coef.slice(i_n8, i_n8 + this.Fr.n8) : this.Fr.zero;

            if (blindingValue !== undefined) {
                b = this.Fr.mul(b, blindingValue);
            }
            if (other) {
                polynomial.coef.set(this.Fr.add(a, b), i_n8);
            } else {
                this.coef.set(this.Fr.add(a, b), i_n8);
            }
        }
        if (other) {
            delete this.coef;
            this.coef = polynomial.coef;
        }
    }

    sub(polynomial, blindingValue) {
        let other = false;

        if (polynomial.length() > this.length()) {
            other = true;
        }

        const thisLength = this.length();
        const polyLength = polynomial.length();
        for (let i = 0; i < Math.max(thisLength, polyLength); i++) {
            const i_n8 = i * this.Fr.n8;

            const a = i < thisLength ? this.coef.slice(i_n8, i_n8 + this.Fr.n8) : this.Fr.zero;
            let b = i < polyLength ? polynomial.coef.slice(i_n8, i_n8 + this.Fr.n8) : this.Fr.zero;

            if (blindingValue !== undefined) {
                b = this.Fr.mul(b, blindingValue);
            }
            if (other) {
                polynomial.coef.set(this.Fr.sub(a, b), i_n8);
            } else {
                this.coef.set(this.Fr.sub(a, b), i_n8);
            }
        }
        if (other) {
            delete this.coef;
            this.coef = polynomial.coef;
        }
    }

    mulScalar(value) {
        for (let i = 0; i < this.length(); i++) {
            const i_n8 = i * this.Fr.n8;

            this.coef.set(this.Fr.mul(this.coef.slice(i_n8, i_n8 + this.Fr.n8), value), i_n8);
        }
    }

    addScalar(value) {
        const currentValue = 0 === this.length() ? this.Fr.zero : this.coef.slice(0, this.Fr.n8);
        this.coef.set(this.Fr.add(currentValue, value), 0);
    }

    subScalar(value) {
        const currentValue = 0 === this.length() ? this.Fr.zero : this.coef.slice(0, this.Fr.n8);
        this.coef.set(this.Fr.sub(currentValue, value), 0);
    }

    // Divide polynomial by X - value
    divByXValue(value) {
        const coefs = new BigBuffer(this.length() * this.Fr.n8);

        coefs.set(this.Fr.zero, (this.length() - 1) * this.Fr.n8);
        coefs.set(this.coef.slice((this.length() - 1) * this.Fr.n8, this.length() * this.Fr.n8), (this.length() - 2) * this.Fr.n8);
        for (let i = this.length() - 3; i >= 0; i--) {
            let i_n8 = i * this.Fr.n8;
            coefs.set(
                this.Fr.add(
                    this.coef.slice(i_n8 + this.Fr.n8, i_n8 + 2 * this.Fr.n8),
                    this.Fr.mul(value, coefs.slice(i_n8 + this.Fr.n8, i_n8 + 2 * this.Fr.n8))
                ),
                i * this.Fr.n8
            );
        }
        if (!this.Fr.eq(
            this.coef.slice(0, this.Fr.n8),
            this.Fr.mul(this.Fr.neg(value), coefs.slice(0, this.Fr.n8))
        )) {
            // throw new Error("Polynomial does not divide");
        }

        this.coef = coefs;
    }

    async divZh() {
        const coefs = new BigBuffer(this.coef.byteLength);

        let domainSize = this.coef.length / 4 / this.Fr.n8;

        if (this.logger) this.logger.debug("dividing T/Z_H");
        for (let i = 0; i < domainSize; i++) {
            const i_n8 = i * this.Fr.n8;
            coefs.set(this.Fr.neg(this.coef.slice(i_n8, i_n8 + this.Fr.n8)), i_n8);
        }

        for (let i = domainSize; i < domainSize * 4; i++) {
            const i_n8 = i * this.Fr.n8;

            const a = this.Fr.sub(
                coefs.slice((i - domainSize) * this.Fr.n8, (i - domainSize) * this.Fr.n8 + this.Fr.n8),
                this.coef.slice(i_n8, i_n8 + this.Fr.n8)
            );
            coefs.set(a, i_n8);
            if (i > (domainSize * 3 - 4)) {
                if (!this.Fr.isZero(a)) {
                    //throw new Error("range_check T Polynomial is not divisible");
                }
            }
        }

        return new Polynomial(coefs, this.Fr, this.logger);
    }

    byX() {
        const coefs = new BigBuffer(this.coef.length + this.Fr.n8);
        coefs.set(this.Fr.zero, 0);
        coefs.set(this.coef, this.Fr.n8);

        this.coef = coefs;
    }

    // Compute a new polynomial f(x^n) from f(x)
    // f(x)   = a_0 + a_1·x + a_2·x^2 + ... + a_j·x^j
    // f(x^n) = a_0 + a_1·x^n + a_2·x^2n + ... + a_j·x^jn
    static async expX(polynomial, n, truncate = false) {
        const Fr = polynomial.Fr;

        if (n < 1) {
            // n == 0 not allowed because it has no sens, but if it's necessary we have to return
            // a zero degree polynomial with a constant coefficient equals to the sum of all the original coefficients
            throw new Error("Compute a new polynomial to a zero or negative number is not allowed");
        } else if (1 === n) {
            return await Polynomial.fromBuffer(polynomial.coef, Fr, polynomial.logger);
        }

        // length is the length of non-constant coefficients
        // if truncate === true, the highest zero coefficients (if exist) will be removed
        const length = truncate ? polynomial.degree() : (polynomial.length() - 1);
        const bufferDst = new BigBuffer((length * n + 1) * Fr.n8);

        // Copy constant coefficient as is because is not related to x
        bufferDst.set(polynomial.coef.slice(0, Fr.n8), 0);

        for (let i = 1; i <= length; i++) {
            const i_sFr = i * Fr.n8;

            const coef = polynomial.coef.slice(i_sFr, i_sFr + Fr.n8);
            bufferDst.set(coef, i_sFr * n);
        }

        return new Polynomial(bufferDst, Fr, polynomial.logger);
    }

    split(numPols, degPols, blindingFactors) {
        if (numPols < 1) {
            throw new Error(`Polynomials can't be split in ${numPols} parts`);
        } else if (1 === numPols) {
            return [this];
        }

        //blinding factors can be void or must have a length of numPols - 1
        if (0 !== blindingFactors.length && blindingFactors.length < numPols - 1) {
            throw new Error(`Blinding factors length must be ${numPols - 1}`);
        }

        const chunkByteLength = (degPols + 1) * this.Fr.n8;
        let res = [];

        // Check polynomial can be split in numChunks parts of chunkSize bytes...
        const numRealPols = Math.ceil((this.degree() + 1) * this.Fr.n8 / chunkByteLength);
        if (numRealPols < numPols) {
            //throw new Error(`Polynomial is short to be split in ${numPols} parts of ${degPols} coefficients each.`);
            for (let i = numRealPols; i < numPols; i++) {
                res[i] = new Polynomial(new Uint8Array(this.Fr.n8), this.Fr, this.logger);
            }
        }

        numPols = Math.min(numPols, numRealPols);
        for (let i = 0; i < numPols; i++) {
            const isLast = (numPols - 1) === i;
            const byteLength = isLast ? this.coef.byteLength - ((numPols - 1) * chunkByteLength) : chunkByteLength + this.Fr.n8;

            res[i] = new Polynomial(new BigBuffer(byteLength), this.Fr, this.logger);
            const fr = i * chunkByteLength;
            const to = isLast ? this.coef.byteLength : (i + 1) * chunkByteLength;
            res[i].coef.set(this.coef.slice(fr, to), 0);

            // Add a blinding factor as higher degree
            if (!isLast) {
                res[i].coef.set(blindingFactors[i], chunkByteLength);
            }

            // Sub blinding factor to the lowest degree
            if (0 !== i) {
                const lowestDegree = this.Fr.sub(res[i].coef.slice(0, this.Fr.n8), blindingFactors[i - 1]);
                res[i].coef.set(lowestDegree, 0);
            }

            if (isLast) {
                res[i].truncate();
            }
        }

        return res;

        // // compute t_low(X)
        // let polTLow = new BigBuffer((chunkSize + 1) * n8r);
        // polTLow.set(t.slice(0, zkey.domainSize * n8r), 0);
        // // Add blinding scalar b_10 as a new coefficient n
        // polTLow.set(ch.b[10], zkey.domainSize * n8r);
        //
        // // compute t_mid(X)
        // let polTMid = new BigBuffer((zkey.domainSize + 1) * n8r);
        // polTMid.set(t.slice(zkey.domainSize * n8r, zkey.domainSize * 2 * n8r), 0);
        // // Subtract blinding scalar b_10 to the lowest coefficient of t_mid
        // const lowestMid = Fr.sub(polTMid.slice(0, n8r), ch.b[10]);
        // polTMid.set(lowestMid, 0);
        // // Add blinding scalar b_11 as a new coefficient n
        // polTMid.set(ch.b[11], zkey.domainSize * n8r);
        //
        // // compute t_high(X)
        // let polTHigh = new BigBuffer((zkey.domainSize + 6) * n8r);
        // polTHigh.set(t.slice(zkey.domainSize * 2 * n8r, (zkey.domainSize * 3 + 6) * n8r), 0);
        // //Subtract blinding scalar b_11 to the lowest coefficient of t_high
        // const lowestHigh = Fr.sub(polTHigh.slice(0, n8r), ch.b[11]);
        // polTHigh.set(lowestHigh, 0);
        //
        // proof.T1 = await expTau(polTLow, "multiexp T1");
        // proof.T2 = await expTau(polTMid, "multiexp T2");
        // proof.T3 = await expTau(polTHigh, "multiexp T3");
    }

    // split2(degPols, blindingFactors) {
    //     let currentDegree = this.degree();
    //     const numFilledPols = Math.ceil((currentDegree + 1) / (degPols + 1));
    //
    //     //blinding factors can be void or must have a length of numPols - 1
    //     if (0 !== blindingFactors.length && blindingFactors.length < numFilledPols - 1) {
    //         throw new Error(`Blinding factors length must be ${numFilledPols - 1}`);
    //     }
    //
    //     const chunkByteLength = (degPols + 1) * this.Fr.n8;
    //
    //     // Check polynomial can be split in numChunks parts of chunkSize bytes...
    //     if (this.coef.byteLength / chunkByteLength <= numFilledPols - 1) {
    //         throw new Error(`Polynomial is short to be split in ${numFilledPols} parts of ${degPols} coefficients each.`);
    //     }
    //
    //     let res = [];
    //     for (let i = 0; i < numFilledPols; i++) {
    //         const isLast = (numFilledPols - 1) === i;
    //         const byteLength = isLast ? (currentDegree + 1) * this.Fr.n8 - ((numFilledPols - 1) * chunkByteLength) : chunkByteLength + this.Fr.n8;
    //
    //         res[i] = new Polynomial(new BigBuffer(byteLength), this.Fr, this.logger);
    //         const fr = i * chunkByteLength;
    //         const to = isLast ? (currentDegree + 1) * this.Fr.n8 : (i + 1) * chunkByteLength;
    //         res[i].coef.set(this.coef.slice(fr, to), 0);
    //
    //         // Add a blinding factor as higher degree
    //         if (!isLast) {
    //             res[i].coef.set(blindingFactors[i], chunkByteLength);
    //         }
    //
    //         // Sub blinding factor to the lowest degree
    //         if (0 !== i) {
    //             const lowestDegree = this.Fr.sub(res[i].coef.slice(0, this.Fr.n8), blindingFactors[i - 1]);
    //             res[i].coef.set(lowestDegree, 0);
    //         }
    //     }
    //
    //     return res;
    // }

    // merge(pols, overlap = true) {
    //     let length = 0;
    //     for (let i = 0; i < pols.length; i++) {
    //         length += pols[i].length();
    //     }
    //
    //     if (overlap) {
    //         length -= pols.length - 1;
    //     }
    //
    //     let res = new Polynomial(new BigBuffer(length * this.Fr.n8));
    //     for (let i = 0; i < pols.length; i++) {
    //         const byteLength = pols[i].coef.byteLength;
    //         if (0 === i) {
    //             res.coef.set(pols[i].coef, 0);
    //         } else {
    //
    //         }
    //     }
    //
    //     return res;
    // }

    truncate() {
        const deg = this.degree();
        if (deg + 1 < this.coef.byteLength / this.Fr.n8) {
            const newCoefs = new BigBuffer((deg + 1) * this.Fr.n8);
            newCoefs.set(this.coef.slice(0, (deg + 1) * this.Fr.n8), 0);
            this.coef = newCoefs;
        }
    }
}