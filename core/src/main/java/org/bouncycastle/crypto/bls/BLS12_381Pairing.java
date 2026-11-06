package org.bouncycastle.crypto.bls;

import java.math.BigInteger;

import org.bouncycastle.math.ec.ECPoint;

/**
 * Optimal ate pairing on BLS12-381: a bilinear, non-degenerate map
 * {@code e: G1 x G2 -> GT} where GT is the order-r subgroup of Fp^12 ^*.
 * <ul>
 *   <li>The Miller loop keeps {@code T} on the twist {@code E'(Fp^2)} in
 *       Jacobian coordinates, so no step inverts, and evaluates each line at
 *       {@code P} directly as a sparse Fp^12 element with three non-zero Fp^2
 *       coefficients, multiplied into the accumulator by
 *       {@link Fp12Element#mulBy014}. The lines are scaled by Fp^2 and Fp^4
 *       factors (the twist's {@code w^3} and the Jacobian denominators); every
 *       such factor lies in a proper subfield of Fp^12, which the final
 *       exponentiation maps to 1, so the pairing value is unchanged.</li>
 *   <li>The final exponentiation splits into the Frobenius-based easy part and
 *       the hard part described at {@code finalExponentiation}.</li>
 * </ul>
 */
public class BLS12_381Pairing
{
    /** BLS12-381 trace parameter |x|: 0xd201000000010000. The actual x is negative. */
    private static final BigInteger ABS_X = new BigInteger("d201000000010000", 16);

    /**
     * Hard-part exponent {@code (p^4 - p^2 + 1) / r} (~1269 bits), applied
     * after the Frobenius-based easy part. The full final exponent
     * {@code (p^12 - 1) / r} (~4317 bits) factors as
     * {@code (p^6 - 1) * (p^2 + 1) * (p^4 - p^2 + 1) / r}; the easy part
     * computes {@code f^((p^6 - 1)(p^2 + 1))} essentially for free using
     * conjugation and one Frobenius&sup2; application, leaving only this
     * shorter exponent, which {@link #hardPart} applies by an addition chain
     * in x (kept for the tests that check that chain against a direct power).
     */
    public static final BigInteger HARD_EXPONENT = Fp2Element.P.pow(4)
        .subtract(Fp2Element.P.pow(2))
        .add(BigInteger.ONE)
        .divide(BLS12_381G1.ORDER);

    /** BLS12-381 trace parameter x = -{@link #ABS_X}. */
    private static final BigInteger X = ABS_X.negate();

    /**
     * {@code (x - 1)^2 / 3}, ~126 bits, used by the cube-factor-free hard
     * part. The division is exact for the BLS12 family by construction.
     */
    private static final BigInteger X_MINUS_1_SQ_OVER_3;

    static
    {
        BigInteger xMinus1Sq = X.subtract(BigInteger.ONE).pow(2);
        BigInteger[] divRem = xMinus1Sq.divideAndRemainder(BigInteger.valueOf(3));
        if (divRem[1].signum() != 0)
        {
            throw new IllegalStateException("(x - 1)^2 must be divisible by 3 for BLS12 family");
        }
        X_MINUS_1_SQ_OVER_3 = divRem[0];
    }

    private BLS12_381Pairing()
    {
    }

    /**
     * Compute the optimal ate pairing {@code e(P, Q)} for {@code P} on
     * BLS12-381 G1 and {@code Q} on BLS12-381 G2.
     *
     * @param g1 a point on the BLS12-381 G1 curve. Must be in the
     *           prime-order subgroup; this method does not subgroup-check.
     * @param g2 a point on the BLS12-381 G2 curve. Must be in the
     *           prime-order subgroup; this method does not subgroup-check.
     * @return {@code e(P, Q)} as an Fp^12 element in the order-r subgroup of
     *         Fp^12 ^*. Returns 1 if either input is the point at infinity.
     */
    public static Fp12Element pair(ECPoint g1, BLS12_381G2Point g2)
    {
        return multiPair(new ECPoint[]{g1}, new BLS12_381G2Point[]{g2});
    }

    /**
     * Multi-pairing: compute the product
     * {@code e(P_0, Q_0) * e(P_1, Q_1) * ... * e(P_{n-1}, Q_{n-1})} with a
     * single shared Miller loop and a single final exponentiation. This
     * cuts a 2-pairing verification (e.g. BLS signature verify) to
     * roughly the cost of one pair() call, since the dominant final
     * exponentiation is performed only once.
     * <p>
     * Pairs whose G1 or G2 component is the point at infinity are skipped
     * (their pairing value is 1, identity in GT). If all pairs are skipped,
     * the result is {@link Fp12Element#ONE}.
     *
     * @param g1Points G1 inputs.
     * @param g2Points G2 inputs; must be the same length as {@code g1Points}.
     * @return the product of pairings as an element of GT.
     * @throws IllegalArgumentException if the arrays differ in length.
     */
    public static Fp12Element multiPair(ECPoint[] g1Points, BLS12_381G2Point[] g2Points)
    {
        if (g1Points.length != g2Points.length)
        {
            throw new IllegalArgumentException("g1 / g2 arrays must be the same length");
        }

        // Skip pairs with a point at infinity: their pairing value is 1.
        int n = 0;
        for (int i = 0; i < g1Points.length; ++i)
        {
            if (!g1Points[i].isInfinity() && !g2Points[i].isInfinity())
            {
                n++;
            }
        }
        if (n == 0)
        {
            return Fp12Element.ONE;
        }

        BigInteger[] negXP = new BigInteger[n];
        BigInteger[] yP = new BigInteger[n];
        Fp2Element[] xQ = new Fp2Element[n];
        Fp2Element[] yQ = new Fp2Element[n];
        Fp2Element[][] t = new Fp2Element[n][];

        int k = 0;
        for (int i = 0; i < g1Points.length; ++i)
        {
            if (g1Points[i].isInfinity() || g2Points[i].isInfinity())
            {
                continue;
            }
            ECPoint normalised = g1Points[i].normalize();
            negXP[k] = normalised.getAffineXCoord().toBigInteger().negate().mod(Fp2Element.P);
            yP[k] = normalised.getAffineYCoord().toBigInteger();
            xQ[k] = g2Points[i].x();
            yQ[k] = g2Points[i].y();
            t[k] = new Fp2Element[]{ xQ[k], yQ[k], Fp2Element.ONE };
            k++;
        }

        Fp12Element f = Fp12Element.ONE;
        int hiBit = ABS_X.bitLength() - 1;
        for (int i = hiBit - 1; i >= 0; --i)
        {
            f = f.square();
            for (int j = 0; j < n; ++j)
            {
                f = doublingStep(f, t[j], negXP[j], yP[j]);
            }
            if (ABS_X.testBit(i))
            {
                for (int j = 0; j < n; ++j)
                {
                    f = additionStep(f, t[j], xQ[j], yQ[j], negXP[j], yP[j]);
                }
            }
        }

        return finalExponentiation(f.conjugate());
    }

    /**
     * Frobenius-based final exponentiation in two stages, producing exactly
     * {@code f^((p^12 - 1) / r)} — the canonical pairing value, with no
     * cube factor.
     * <p>
     * <b>Easy part</b> {@code f^((p^6 - 1)(p^2 + 1))}: computed as
     * {@code conjugate(f) * inverse(f)} (one Fp&sup12; inverse) followed by
     * {@code result * frobeniusSquared(result)} (one Frobenius&sup2; +
     * one mul). After the easy part, the result lives in the cyclotomic
     * subgroup of Fp&sup12;, where {@code conjugate} equals inversion.
     * <p>
     * <b>Hard part</b> {@code easy^((p^4 - p^2 + 1) / r)} via the
     * Hayashida-Hayasaka-Teruya 2020 decomposition for the BLS12 family
     * with k = 12, rearranged to factor out the +3 term that would
     * otherwise leave a cube factor:
     * <pre>
     *   3 * Phi_12(p) / r = (x - 1)^2 * (x + p) * (x^2 + p^2 - 1) + 3
     *   =>  hard          = [(x - 1)^2 / 3] * (x + p) * (x^2 + p^2 - 1) + 1
     * </pre>
     * For the BLS12 family, {@code (x - 1)^2} is divisible by 3 by
     * construction (the curve parameterisation requires it for
     * {@code p(x)} to have integer coefficients), so the {@code /3}
     * collapses into a precomputable ~126-bit integer that is applied
     * via a single cyclotomic power by {@link #X_MINUS_1_SQ_OVER_3}.
     * The remaining structure is the same as before, but the final
     * {@code +3} becomes {@code +1}, so we multiply by {@code f} rather
     * than {@code f^3}. Net cost is slightly higher than the cube-factored
     * variant (a ~126-bit power replaces three exp-by-|x| ops + a square)
     * but produces a pairing value byte-comparable against any other
     * BLS12-381 implementation (blst, mcl, zkcrypto, …).
     */
    private static Fp12Element finalExponentiation(Fp12Element f)
    {
        Fp12Element f1 = f.conjugate().mul(f.inverse());            // f^(p^6 - 1)
        Fp12Element easy = f1.frobeniusSquared().mul(f1);            // f^((p^6 - 1)(p^2 + 1))
        return hardPart(easy);
    }

    /**
     * The hard part of the final exponentiation, exposed for cross-package
     * layered testing (the test classes live in
     * {@code org.bouncycastle.crypto.hash2curve.test} and need direct access
     * to the easy/hard split for KAT comparison against reference outputs).
     * Not part of the intended public API of this class — production callers
     * should use {@link #pair} / {@link #multiPair}.
     * <p>
     * {@code f} must lie in the cyclotomic subgroup, as the output of the easy
     * part does: every power here squares with
     * {@code Fp12Element.cyclotomicSquare}, which is wrong outside it.
     */
    public static Fp12Element hardPart(Fp12Element f)
    {
        // g0 = f^((x - 1)^2 / 3), via a direct ~126-bit power. This is
        // the only step that differs from the 3*hard variant; it replaces
        // f^((x - 1)^2) and absorbs the /3.
        Fp12Element g0 = f.cyclotomicPow(X_MINUS_1_SQ_OVER_3);

        // g1 = g0^(x + p) = g0^x * g0^p.
        Fp12Element g1 = expByX(g0).mul(g0.frobenius());

        // g2 = g1^(x^2 + p^2 - 1) = g1^(x^2) * g1^(p^2) * g1^(-1).
        Fp12Element g1XSq = expByX(expByX(g1));
        Fp12Element g2 = g1XSq.mul(g1.frobeniusSquared()).mul(g1.conjugate());

        // result = g2 * f (the +1 term, vs +3 -> f^3 in the cube-factor variant).
        return g2.mul(f);
    }

    /**
     * In the cyclotomic subgroup of Fp&sup12;, inversion is conjugation. So
     * {@code f^x = f^(-|x|) = conjugate(f^|x|)} for our negative parameter x.
     */
    private static Fp12Element expByX(Fp12Element f)
    {
        return f.cyclotomicPow(ABS_X).conjugate();
    }

    /**
     * Doubles {@code T = (X, Y, Z)} in place (Jacobian, dbl-2009-l) and returns
     * {@code f} times the tangent line at {@code T} evaluated at
     * {@code P = (xP, yP)}.
     * <p>
     * On the twist the tangent has slope {@code 3x^2 / 2y}; after untwisting and
     * multiplying by {@code w^3} its value at {@code P} is
     * {@code (lambda*x - y) - lambda*xP*v + yP*v*w}. Scaled by
     * {@code 2*Y*Z^3}, with {@code x = X/Z^2}, {@code y = Y/Z^3}, the three
     * coefficients are {@code 3X^3 - 2Y^2}, {@code -3X^2*Z^2*xP} and
     * {@code 2Y*Z^3*yP}.
     */
    private static Fp12Element doublingStep(Fp12Element f, Fp2Element[] t, BigInteger negXP, BigInteger yP)
    {
        Fp2Element x1 = t[0], y1 = t[1], z1 = t[2];
        if (y1.isZero())
        {
            throw new ArithmeticException("Miller loop doubled a point of order 2: G2 input not in the subgroup");
        }
        Fp2Element a = x1.square();
        Fp2Element b = y1.square();
        Fp2Element c = b.square();
        Fp2Element z1z1 = z1.square();
        Fp2Element d = x1.add(b).square().sub(a).sub(c);
        d = d.add(d);
        Fp2Element e = a.add(a).add(a);
        Fp2Element x3 = e.square().sub(d.add(d));
        Fp2Element c8 = c.add(c);
        c8 = c8.add(c8);
        c8 = c8.add(c8);
        Fp2Element y3 = e.mul(d.sub(x3)).sub(c8);
        Fp2Element yz = y1.mul(z1);
        Fp2Element z3 = yz.add(yz);

        Fp2Element l0 = x1.mul(e).sub(b.add(b));
        Fp2Element l1 = e.mul(z1z1).mulFp(negXP);
        Fp2Element l4 = z3.mul(z1z1).mulFp(yP);

        t[0] = x3;
        t[1] = y3;
        t[2] = z3;
        return f.mulBy014(l0, l1, l4);
    }

    /**
     * Adds the affine {@code Q = (xQ, yQ)} to {@code T = (X, Y, Z)} in place
     * (mixed Jacobian addition, madd-2007-bl) and returns {@code f} times the
     * line through {@code T} and {@code Q} evaluated at {@code P}.
     * <p>
     * With {@code lambda = (yQ - y) / (xQ - x)} the untwisted line, multiplied
     * by {@code w^3}, is {@code (lambda*xQ - yQ) - lambda*xP*v + yP*v*w}; scaled
     * by {@code Z3 = 2*H*Z} (so {@code lambda * Z3 = r}) its coefficients are
     * {@code r*xQ - yQ*Z3}, {@code -r*xP} and {@code Z3*yP}.
     */
    private static Fp12Element additionStep(Fp12Element f, Fp2Element[] t, Fp2Element xQ, Fp2Element yQ,
        BigInteger negXP, BigInteger yP)
    {
        Fp2Element x1 = t[0], y1 = t[1], z1 = t[2];
        Fp2Element z1z1 = z1.square();
        Fp2Element u2 = xQ.mul(z1z1);
        Fp2Element s2 = yQ.mul(z1).mul(z1z1);
        Fp2Element h = u2.sub(x1);
        if (h.isZero())
        {
            // T = +-Q: impossible for Q of prime order r, since T = [k]Q with 1 < k < r.
            throw new ArithmeticException("Miller loop reached T = +-Q: G2 input not in the subgroup");
        }
        Fp2Element hh = h.square();
        Fp2Element i4 = hh.add(hh);
        i4 = i4.add(i4);
        Fp2Element j = h.mul(i4);
        Fp2Element r = s2.sub(y1);
        r = r.add(r);
        Fp2Element v = x1.mul(i4);
        Fp2Element x3 = r.square().sub(j).sub(v.add(v));
        Fp2Element y1j = y1.mul(j);
        Fp2Element y3 = r.mul(v.sub(x3)).sub(y1j.add(y1j));
        Fp2Element z3 = z1.add(h).square().sub(z1z1).sub(hh);

        Fp2Element l0 = r.mul(xQ).sub(yQ.mul(z3));
        Fp2Element l1 = r.mulFp(negXP);
        Fp2Element l4 = z3.mulFp(yP);

        t[0] = x3;
        t[1] = y3;
        t[2] = z3;
        return f.mulBy014(l0, l1, l4);
    }
}
