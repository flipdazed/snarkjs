// SPDX-License-Identifier: GPL-3.0
/*
    Copyright 2021 0KIMS association.

    This file is generated with [snarkJS](https://github.com/iden3/snarkjs).

    snarkJS is a free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    snarkJS is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
    License for more details.

    You should have received a copy of the GNU General Public License
    along with snarkJS. If not, see <https://www.gnu.org/licenses/>.
*/


pragma solidity >=0.7.0 <0.9.0;

contract PlonkVerifier {

    uint32 constant n         = 32;
    uint16 constant nPublic   = 1;
    uint16 constant nLagrange = 1;

    // Verification Key data
    uint256 constant k1   = 2;
    uint256 constant k2   = 3;
    uint256 constant w1   = 4419234939496763621076330863786513495701855246241724391626358375488475697872;
    uint256 constant w3   = 21888242871839275217838484774961031246154997185409878258781734729429964517155;
    uint256 constant w3_2 = 4407920970296243842393367215006156084916469457145843978461;
    uint256 constant w4   = 21888242871839275217838484774961031246007050428528088939761107053157389710902;
    uint256 constant w4_2 = 21888242871839275222246405745257275088548364400416034343698204186575808495616;
    uint256 constant w4_3 = 4407920970296243842541313971887945403937097133418418784715;
    uint256 constant wr   = 9222527969605388450625148037496647087331675164191659244434925070698893435503;
    uint256 constant X2x1 = 18029695676650738226693292988307914797657423701064905010927197838374790804409;
    uint256 constant X2x2 = 14583779054894525174450323658765874724019480979794335525732096752006891875705;
    uint256 constant X2y1 = 2140229616977736810657479771656733941598412651537078903776637920509952744750;
    uint256 constant X2y2 = 11474861747383700316476719153975578001603231366361248090558603872215261634898;

    uint256 constant q    = 21888242871839275222246405745257275088548364400416034343698204186575808495617;
    uint256 constant qf   = 21888242871839275222246405745257275088696311157297823662689037894645226208583;

    uint256 constant G1x  = 1;
    uint256 constant G1y  = 2;
    uint256 constant G2x1 = 10857046999023057135944570762232829481370756359578518086990519993285655852781;
    uint256 constant G2x2 = 11559732032986387107991004021392285783925812861821192530917403151452391805634;
    uint256 constant G2y1 = 8495653923123431417604973247489272438418190587263600148770280649306958101930;
    uint256 constant G2y2 = 4082367875863433681332203403145435568316851327593401208105741076214120093531;

    // Proof data
    uint16 constant pC1 = 32;
    uint16 constant pC2 = 96;
    uint16 constant pW1 = 160;
    uint16 constant pW2 = 224;
    uint16 constant pEval_ql  = 288;
    uint16 constant pEval_qr  = 320;
    uint16 constant pEval_qm  = 352;
    uint16 constant pEval_qo  = 384;
    uint16 constant pEval_qc  = 416;
    uint16 constant pEval_s1  = 448;
    uint16 constant pEval_s2  = 480;
    uint16 constant pEval_s3  = 512;
    uint16 constant pEval_a   = 544;
    uint16 constant pEval_b   = 576;
    uint16 constant pEval_c   = 608;
    uint16 constant pEval_z   = 640;
    uint16 constant pEval_zw  = 672;
    uint16 constant pEval_t1w = 704;
    uint16 constant pEval_t2w = 736;
    uint16 constant pEval_inv = 768;

    // Memory data
    // Challenges
    uint16 constant pAlpha  = 0;
    uint16 constant pBeta   = 32;
    uint16 constant pGamma  = 64;
    uint16 constant pY      = 96;
    uint16 constant pXiSeed = 128;
    uint16 constant pXiSeed2= 160;
    uint16 constant pXi     = 192;

    // Roots
    uint16 constant pH1w4_0 = 224;
    uint16 constant pH1w4_1 = 256;
    uint16 constant pH1w4_2 = 288;
    uint16 constant pH1w4_3 = 320;

    uint16 constant pH2w3_0 = 352;
    uint16 constant pH2w3_1 = 384;
    uint16 constant pH2w3_2 = 416;

    uint16 constant pH3w3_0 = 448;
    uint16 constant pH3w3_1 = 480;
    uint16 constant pH3w3_2 = 512;

    uint16 constant pBetaXi = 544;

    uint16 constant pPi     = 576;
    uint16 constant pTmp    = 608;
    uint16 constant pQuo    = 640;
    uint16 constant pR1     = 672;
    uint16 constant pR2     = 704;

    uint16 constant pF      = 736;  // 64 bytes
    uint16 constant pE      = 800;  // 64 bytes
    uint16 constant pJ      = 864;  // 64 bytes
    uint16 constant pA      = 928; // 64 bytes

    uint16 constant pZh     = 992;
    // From this point we write all the variables that must compute the inverse using the Montgomery batch inversion
    uint16 constant pZhInv  = 1024;
    uint16 constant pDen    = 1056;
    uint16 constant pRInv   = 1088; // Reserve 10 * 32 bytes to compute R1 and R2 inversions
    
    uint16 constant pEval_l1 = 1408;
    
    
    uint16 constant lastMem = 1440;

    function verifyProof(bytes memory proof, uint[] memory pubSignals) public view returns (bool) {
        assembly {
            ///////
            // Computes the inverse of an array of values
            // See https://vitalik.ca/general/2018/07/21/starks_part_3.html in section where explain fields operations
            //////
            function inverseArray(pProof, pVals, n) {

                let pAux := mload(0x40)     // Point to the next free position
                let pIn := pVals
                let lastPIn := add(pVals, mul(n, 32))  // Read n elemnts
                let acc := mload(pIn)       // Read the first element
                pIn := add(pIn, 32)         // Point to the second element

                for { } lt(pIn, lastPIn) {
                    pAux := add(pAux, 32)
                    pIn := add(pIn, 32)
                }
                {
                    mstore(pAux, acc)
                    acc := mulmod(acc, mload(pIn), q)
                }

                let inv := mload(add(pProof, pEval_inv))
                // Check inv * 1/inv == 1
                if iszero(eq(1, mulmod(acc, inv, q) )) {
                    mstore(0, 0)
                    return(0,0x20)
                }

                acc := inv

                // At this point pAux pint to the next free position we substract 1 to point to the last used
                pAux := sub(pAux, 32)
                // pIn points to the n+1 element, we substract to point to n
                pIn := sub(pIn, 32)
                lastPIn := pVals  // We don't process the first element
                for { } gt(pIn, lastPIn) {
                    pAux := sub(pAux, 32)
                    pIn := sub(pIn, 32)
                }
                {
                    inv := mulmod(acc, mload(pAux), q)
                    acc := mulmod(acc, mload(pIn), q)
                    mstore(pIn, inv)
                }
                // pIn points to first element, we just set it.
                mstore(pIn, acc)
            }

            function checkField(v) {
                if iszero(lt(v, q)) {
                    mstore(0, 0)
                    return(0,0x20)
                }
            }

            function checkInput(pProof) {
                if iszero(eq(mload(pProof), 768 )) {
                    mstore(0, 0)
                    return(0,0x20)
                }

                checkField(mload(add(pProof, pEval_ql)))
                checkField(mload(add(pProof, pEval_qr)))
                checkField(mload(add(pProof, pEval_qm)))
                checkField(mload(add(pProof, pEval_qo)))
                checkField(mload(add(pProof, pEval_qc)))
                checkField(mload(add(pProof, pEval_s1)))
                checkField(mload(add(pProof, pEval_s2)))
                checkField(mload(add(pProof, pEval_s3)))
                checkField(mload(add(pProof, pEval_a)))
                checkField(mload(add(pProof, pEval_b)))
                checkField(mload(add(pProof, pEval_c)))
                checkField(mload(add(pProof, pEval_z)))
                checkField(mload(add(pProof, pEval_zw)))
                checkField(mload(add(pProof, pEval_t1w)))
                checkField(mload(add(pProof, pEval_t2w)))
                checkField(mload(add(pProof, pEval_inv)))

                // Points are checked in the point operations precompiled smart contracts
            }

            function computeChallenges(pProof, pMem, pPublic) {
                // Compute challenge.beta & challenge.gamma
                let challenge
                
                mstore( add(pMem, 1440 ), mload( add( pPublic, 32)))
                
                mstore( add(pMem, 1472 ),  mload( add( pProof, pC1)))
                mstore( add(pMem, 1504 ),  mload( add( pProof, add(pC1, 32))))

                challenge := mod(keccak256(add(pMem, lastMem), 96), q)

                mstore( add(pMem, pBeta), challenge)
                mstore( add(pMem, pGamma), mod(keccak256(add(pMem, pBeta), 32), q))

                // Get xiSeed & xiSeed2
                mstore( add(pMem, lastMem ),  mload( add( pProof, pC2)))
                mstore( add(pMem, 1472 ),  mload( add( pProof, add(pC2, 32))))
                challenge := mod(keccak256(add(pMem, lastMem), 64), q)

                mstore( add(pMem, pXiSeed), challenge)
                mstore( add(pMem, pXiSeed2), mulmod(challenge, challenge, q))

                // Compute roots.S1.h1w4
                mstore( add(pMem, pH1w4_0), mulmod(mload(add(pMem, pXiSeed2)), mload(add(pMem, pXiSeed)), q))
                mstore( add(pMem, pH1w4_1), mulmod(mload(add(pMem, pH1w4_0)), w4, q))
                mstore( add(pMem, pH1w4_2), mulmod(mload(add(pMem, pH1w4_0)), w4_2, q))
                mstore( add(pMem, pH1w4_3), mulmod(mload(add(pMem, pH1w4_0)), w4_3, q))

                // Compute roots.S2.h2w3
                mstore( add(pMem, pH2w3_0), mulmod(mload(add(pMem, pXiSeed2)), mload(add(pMem, pXiSeed2)), q))
                mstore( add(pMem, pH2w3_1), mulmod(mload(add(pMem, pH2w3_0)), w3, q))
                mstore( add(pMem, pH2w3_2), mulmod(mload(add(pMem, pH2w3_0)), w3_2, q))

                // Compute roots.S2.h2w3
                mstore( add(pMem, pH3w3_0), mulmod(mload(add(pMem, pH2w3_0)), wr, q))
                mstore( add(pMem, pH3w3_1), mulmod(mload(add(pMem, pH3w3_0)), w3, q))
                mstore( add(pMem, pH3w3_2), mulmod(mload(add(pMem, pH3w3_0)), w3_2, q))

                let xin := mulmod(mulmod(mload(add(pMem, pH2w3_0)), mload(add(pMem, pH2w3_0)), q), mload(add(pMem, pH2w3_0)), q)
                mstore( add(pMem, pXi), xin)

                // Compute beta*xi and store it in a local variable to use it later
                mstore( add(pMem, pBetaXi), mulmod(mload(add(pMem, pBeta)), mload(add(pMem, pXi)), q))

                // Compute xi^n
                
                    xin:= mulmod(xin, xin, q)
                
                    xin:= mulmod(xin, xin, q)
                
                    xin:= mulmod(xin, xin, q)
                
                    xin:= mulmod(xin, xin, q)
                
                    xin:= mulmod(xin, xin, q)
                

                xin:= mod(add(sub(xin, 1), q), q)
                mstore( add(pMem, pZh), xin)
                mstore( add(pMem, pZhInv), xin)  // We will invert later together with lagrange pols

                // Compute challenge.alpha
                mstore( add(pMem, pAlpha), mod(keccak256(add(pProof, pEval_ql), 480), q))

                // Compute challenge.y
                mstore( add(pMem, pY), mod(keccak256(add(pProof, pW1), 64), q))
            }

            function precomputeF(pMem) {
                let den

                den := addmod(mload(add(pMem, pY)), mod(sub(q, mload(add(pMem, pH2w3_0))), q), q)
                den := mulmod(den, addmod(mload(add(pMem, pY)), mod(sub(q, mload(add(pMem, pH2w3_1))), q), q), q)
                den := mulmod(den, addmod(mload(add(pMem, pY)), mod(sub(q, mload(add(pMem, pH2w3_2))), q), q), q)
                den := mulmod(den, addmod(mload(add(pMem, pY)), mod(sub(q, mload(add(pMem, pH3w3_0))), q), q), q)
                den := mulmod(den, addmod(mload(add(pMem, pY)), mod(sub(q, mload(add(pMem, pH3w3_1))), q), q), q)
                den := mulmod(den, addmod(mload(add(pMem, pY)), mod(sub(q, mload(add(pMem, pH3w3_2))), q), q), q)

                mstore(add(pMem, pDen), den)
            }

            function calcLagrangeItem(pMem, i, n, pX, pRoot) -> result {
                let idx := i
                let max := add(n, 1)
                result := 1
                for { let j := 0 } lt(j, n) { j := add(j, 1) }
                {
                    idx := mod(add(idx, 1), max)

                    result := mulmod(result, addmod(mload(add(pMem, pX)), mod(sub(q, mload(add(pMem, add(pRoot, mul(idx, 32))))), q), q), q)
                }
            }

            function precomputeR1(pMem) {
                for { let i := 0 } lt(i, 4) { i := add(i, 1) }
                {
                    mstore(add(pMem, add(pRInv, mul(i, 32))), calcLagrangeItem(pMem, i, 3, add(pH1w4_0, mul(i, 32)), pH1w4_0))
                }
            }

            function precomputeR2(pMem) {
                for { let i := 0 } lt(i, 6) { i := add(i, 1) }
                {
                    mstore(add(pMem, add(pRInv, add(128, mul(i, 32)))), calcLagrangeItem(pMem, i, 5, add(pH2w3_0, mul(i, 32)), pH2w3_0))
                }
            }

            function precomputeLagrange(pMem) {
                let w := 1
                
                mstore(add(pMem, pEval_l1), mulmod(n, mod(add(sub(mload(add(pMem, pXi)), w), q), q), q))
                
                
            }

            function computeInversions(pProof, pMem) {
                inverseArray(pProof, add(pMem, pZhInv), 13)
            }

            function computeLagrange(pMem) {
                let zh := mload(add(pMem, pZh))
                let w := 1
                
                
                mstore(add(pMem, pEval_l1 ), mulmod(mload(add(pMem, pEval_l1 )), zh, q))
                
                
                
            }

            function computePi(pMem, pPub) {
                let pi := 0

                
                pi := mod(
                    add(
                        sub(
                            pi,
                            mulmod(
                                mload(add(pMem, pEval_l1)),
                                mload(add(pPub, 32)),
                                q
                            )
                        ),
                        q
                    ),
                    q
                )
                

                mstore(add(pMem, pPi), pi)
            }

            function computeR1(pProof, pMem) {
                let t0

                t0 := mulmod(mload(add(pProof, pEval_ql)), mload(add(pProof, pEval_a)), q)
                t0 := addmod(t0, mulmod(mload(add(pProof, pEval_qr)), mload(add(pProof, pEval_b)), q) ,q)
                t0 := addmod(t0, mulmod(mload(add(pProof, pEval_qm)),
                mulmod(mload(add(pProof, pEval_a)), mload(add(pProof, pEval_b)), q)
                , q) ,q)
                t0 := addmod(t0, mulmod(mload(add(pProof, pEval_qo)), mload(add(pProof, pEval_c)), q) ,q)
                t0 := addmod(t0, mload(add(pProof, pEval_qc)) ,q)
                t0 := addmod(t0, mload(add(pMem, pPi)), q)
                t0 := mulmod(t0, mload(add(pMem, pZhInv)), q)

                let res
                for { let i := 0 } lt(i, 4) { i := add(i, 1) }
                {
                    let c1Value := mload(add(pProof, pEval_a))
                    let h1w4 := mload(add(pMem, add(pH1w4_0, mul(i, 32))))

                    c1Value := addmod(c1Value, mulmod(h1w4, mload(add(pProof, pEval_b)), q), q)
                    let square := mulmod(h1w4, h1w4, q)
                    c1Value := addmod(c1Value, mulmod(square, mload(add(pProof, pEval_c)), q), q)
                    c1Value := addmod(c1Value, mulmod(mulmod(square, h1w4, q), t0, q), q)

                    let lagrange := calcLagrangeItem(pMem, i, 3, pY, pH1w4_0)
                    lagrange := mulmod(lagrange, mload(add(pMem, add(pRInv, mul(i, 32)))), q)

                    res := addmod(res, mulmod(c1Value, lagrange, q), q)
                }

                mstore(add(pMem, pR1), res)
            }

            function computeR2(pProof, pMem) {
                let t2
                t2 := addmod(mload(add(pProof, pEval_a)), addmod(mload(add(pMem, pBetaXi)), mload(add(pMem, pGamma)), q) ,q)
                t2 := mulmod(t2,
                            addmod(mload(add(pProof, pEval_b)),
                            addmod(mulmod(mload(add(pMem, pBetaXi)), k1, q), mload(add(pMem, pGamma)), q) ,q), q)
                t2 := mulmod(t2,
                            addmod(mload(add(pProof, pEval_c)),
                            addmod(mulmod(mload(add(pMem, pBetaXi)), k2, q), mload(add(pMem, pGamma)), q) ,q), q)
                t2 := mulmod(t2, mload(add(pProof, pEval_z)), q)

                let t1 //Let's use this variable as a temporal variable to save one local
                t1 := addmod(mload(add(pProof, pEval_a)), addmod(mulmod(mload(add(pMem, pBeta)), mload(add(pProof, pEval_s1)), q), mload(add(pMem, pGamma)), q) ,q)
                t1 := mulmod(t1,
                      addmod(mload(add(pProof, pEval_b)), addmod(mulmod(mload(add(pMem, pBeta)), mload(add(pProof, pEval_s2)), q), mload(add(pMem, pGamma)), q) ,q), q)
                t1 := mulmod(t1,
                      addmod(mload(add(pProof, pEval_c)), addmod(mulmod(mload(add(pMem, pBeta)), mload(add(pProof, pEval_s3)), q), mload(add(pMem, pGamma)), q) ,q), q)
                t1 := mulmod(t1, mload(add(pProof, pEval_zw)), q)

                t2:= addmod(t2, mod(sub(q, t1), q), q)
                t2 := mulmod(t2, mload(add(pMem, pZhInv)), q)

                // Compute T1(xi)
                t1 := sub(mload(add(pProof, pEval_z)), 1)
                t1 := mulmod(t1, mload(add(pMem, pEval_l1)) ,q)
                t1 := mulmod(t1, mload(add(pMem, pZhInv)) ,q)

                let res
                for { let i := 0 } lt(i, 6) { i := add(i, 1) }
                {
                    let hw := mload(add(pMem, add(pH2w3_0, mul(i, 32))))

                    let c2Value
                    if lt(i, 3) {
                        c2Value := addmod(mload(add(pProof, pEval_z)), mulmod(hw, t1, q), q)
                        c2Value := addmod(c2Value, mulmod(mulmod(hw, hw, q), t2, q), q)
                    }
                    if gt(i, 2) {
                        c2Value := addmod(mload(add(pProof, pEval_zw)), mulmod(hw, mload(add(pProof, pEval_t1w)), q), q)
                        c2Value := addmod(c2Value, mulmod(mulmod(hw, hw, q), mload(add(pProof, pEval_t2w)), q), q)
                    }

                    let lagrange := calcLagrangeItem(pMem, i, 5, pY, pH2w3_0)
                    lagrange := mulmod(lagrange, mload(add(pMem, add(pRInv,add(128, mul(i, 32))))), q)

                    res := addmod(res, mulmod(c2Value, lagrange, q), q)
                }

                mstore(add(pMem, pR2), res)
            }

            function g1_set(pR, pP) {
                mstore(pR, mload(pP))
                mstore(add(pR, 32), mload(add(pP,32)))
            }

            function g1_acc(pR, pP) {
                let mIn := mload(0x40)
                mstore(mIn, mload(pR))
                mstore(add(mIn,32), mload(add(pR, 32)))
                mstore(add(mIn,64), mload(pP))
                mstore(add(mIn,96), mload(add(pP, 32)))

                let success := staticcall(sub(gas(), 2000), 6, mIn, 128, pR, 64)

                if iszero(success) {
                    mstore(0, 0)
                    return(0,0x20)
                }
            }

            function g1_mulAcc(pR, pP, s) {
                let success
                let mIn := mload(0x40)
                mstore(mIn, mload(pP))
                mstore(add(mIn,32), mload(add(pP, 32)))
                mstore(add(mIn,64), s)

                success := staticcall(sub(gas(), 2000), 7, mIn, 96, mIn, 64)

                if iszero(success) {
                    mstore(0, 0)
                    return(0,0x20)
                }

                mstore(add(mIn,64), mload(pR))
                mstore(add(mIn,96), mload(add(pR, 32)))

                success := staticcall(sub(gas(), 2000), 6, mIn, 128, pR, 64)

                if iszero(success) {
                    mstore(0, 0)
                    return(0,0x20)
                }
            }

            function g1_mulAccC(pR, x, y, s) {
                let success
                let mIn := mload(0x40)
                mstore(mIn, x)
                mstore(add(mIn,32), y)
                mstore(add(mIn,64), s)

                success := staticcall(sub(gas(), 2000), 7, mIn, 96, mIn, 64)

                if iszero(success) {
                    mstore(0, 0)
                    return(0,0x20)
                }

                mstore(add(mIn,64), mload(pR))
                mstore(add(mIn,96), mload(add(pR, 32)))

                success := staticcall(sub(gas(), 2000), 6, mIn, 128, pR, 64)

                if iszero(success) {
                    mstore(0, 0)
                    return(0,0x20)
                }
            }

            function g1_mulSetC(pR, x, y, s) {
                let success
                let mIn := mload(0x40)
                mstore(mIn, x)
                mstore(add(mIn,32), y)
                mstore(add(mIn,64), s)

                success := staticcall(sub(gas(), 2000), 7, mIn, 96, pR, 64)

                if iszero(success) {
                    mstore(0, 0)
                    return(0,0x20)
                }
            }

            function computeF(pProof, pMem) {
                let num

                num := addmod(mload(add(pMem, pY)), mod(sub(q, mload(add(pMem, pH1w4_0))), q), q)
                num := mulmod(num , addmod(mload(add(pMem, pY)), mod(sub(q, mload(add(pMem, pH1w4_1))), q), q), q)
                num := mulmod(num , addmod(mload(add(pMem, pY)), mod(sub(q, mload(add(pMem, pH1w4_2))), q), q), q)
                num := mulmod(num , addmod(mload(add(pMem, pY)), mod(sub(q, mload(add(pMem, pH1w4_3))), q), q), q)

                mstore(add(pMem, pTmp), num)

                mstore(add(pMem, pQuo), mulmod(mload(add(pMem, pAlpha)), mulmod(num, mload(add(pMem, pDen)), q), q))

                let p := add(pMem, pF)
                g1_mulAcc(p, add(pProof, pC2), mload(add(pMem, pQuo)))
                g1_acc(p, add(pProof, pC1))
            }

            function computeE(pProof, pMem) {
                let e := addmod(mload(add(pMem, pR1)), mulmod(mload(add(pMem, pQuo)), mload(add(pMem, pR2)),q) , q)
                g1_mulAccC(add(pMem, pE), G1x, G1y, e)
            }

            function computeJ(pProof, pMem) {
                g1_mulAcc(add(pMem, pJ), add(pProof, pW1), mload(add(pMem, pTmp)))
            }

            function checkPairing(pProof, pMem) -> isOk {
                let mIn := mload(0x40)

                // First pairing value
                // Compute -E
                mstore(add(add(pMem, pE), 32), mod(sub(qf, mload(add(add(pMem, pE), 32))), qf))
                // Compute -J
                mstore(add(add(pMem, pJ), 32), mod(sub(qf, mload(add(add(pMem, pJ), 32))), qf))
                // A = F - E - J + y·W2
                g1_acc(add(pMem, pA), add(pMem, pF))
                g1_acc(add(pMem, pA), add(pMem, pE))
                g1_acc(add(pMem, pA), add(pMem, pJ))
                g1_mulAcc(add(pMem, pA), add(pProof, pW2), mload(add(pMem, pY)))

                mstore(mIn, mload(add(pMem, pA)))
                mstore(add(mIn,32), mload(add(add(pMem, pA), 32)))

                // Second pairing value
                mstore(add(mIn,64), G2x2)
                mstore(add(mIn,96), G2x1)
                mstore(add(mIn,128), G2y2)
                mstore(add(mIn,160), G2y1)

                // Third pairing value
                // Compute -W2
                mstore(add(mIn, 192), mload(add(pProof, pW2)))
                let s := mload(add(add(pProof, pW2), 32))
                s := mod(sub(qf, s), qf)
                mstore(add(mIn,224), s)

                // Fourth pairing value
                mstore(add(mIn,256), X2x2)
                mstore(add(mIn,288), X2x1)
                mstore(add(mIn,320), X2y2)
                mstore(add(mIn,352), X2y1)

                let success := staticcall(sub(gas(), 2000), 8, mIn, 384, mIn, 0x20)

                isOk := and(success, mload(mIn))
            }

            let pMem := mload(0x40)
            mstore(0x40, add(pMem, lastMem))

            checkInput(proof)
            computeChallenges(proof, pMem, pubSignals)

            precomputeF(pMem)
            precomputeR1(pMem)
            precomputeR2(pMem)
            precomputeLagrange(pMem)

            computeInversions(proof, pMem)

            computeLagrange(pMem)
            computePi(pMem, pubSignals)
            computeR1(proof, pMem)
            computeR2(proof, pMem)

            computeF(proof, pMem)
            computeE(proof, pMem)
            computeJ(proof, pMem)
            let isValid := checkPairing(proof, pMem)

            mstore(0x40, sub(pMem, lastMem))
            mstore(0, isValid)
            return(0,0x20)
        }

    }
}
