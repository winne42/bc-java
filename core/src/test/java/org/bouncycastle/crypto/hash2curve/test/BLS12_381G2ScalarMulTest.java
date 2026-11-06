package org.bouncycastle.crypto.hash2curve.test;

import java.math.BigInteger;
import java.security.SecureRandom;

import junit.framework.TestCase;
import org.bouncycastle.crypto.bls.BLS12_381G1;
import org.bouncycastle.crypto.bls.BLS12_381G2;
import org.bouncycastle.crypto.bls.BLS12_381G2HashToCurve;
import org.bouncycastle.crypto.bls.BLS12_381G2Point;
import org.bouncycastle.crypto.bls.Fp2Element;

/**
 * Checks both G2 scalar multiplications, which run in Jacobian coordinates,
 * against a plain affine double-and-add built from the public {@code add} and
 * {@code doublePoint}. Covers points in and outside the prime-order subgroup,
 * negative scalars, and the scalars that drive the mixed addition into its
 * special cases: {@code r} (the accumulator reaches {@code -P}, so the final
 * addition gives infinity) and {@code r + 2} (it reaches {@code P}, so the final
 * addition is a doubling).
 */
public class BLS12_381G2ScalarMulTest
    extends TestCase
{
    private static final BigInteger R = BLS12_381G1.ORDER;
    private static final SecureRandom RNG = new SecureRandom(new byte[]{ 11 });

    private static BLS12_381G2Point affineMultiply(BLS12_381G2Point p, BigInteger k)
    {
        if (k.signum() < 0)
        {
            return affineMultiply(p.negate(), k.negate());
        }
        BLS12_381G2Point result = BLS12_381G2Point.INFINITY;
        BLS12_381G2Point addend = p;
        for (int i = 0; i < k.bitLength(); ++i)
        {
            if (k.testBit(i))
            {
                result = result.add(addend);
            }
            addend = addend.doublePoint();
        }
        return result;
    }

    private static BLS12_381G2Point nonSubgroupPoint()
    {
        BLS12_381G2HashToCurve h = new BLS12_381G2HashToCurve(new byte[]{ 'T', 'E', 'S', 'T' });
        return h.mapToCurveAndIso3(Fp2Element.of(new BigInteger(380, RNG), new BigInteger(380, RNG)));
    }

    private static void check(BLS12_381G2Point p, BigInteger k, boolean constantTime)
    {
        BLS12_381G2Point expected = affineMultiply(p, k);
        assertEquals("multiply by " + k, expected, p.multiply(k));
        if (constantTime)
        {
            assertEquals("constantTimeMultiply by " + k, expected, p.constantTimeMultiply(k));
        }
    }

    public void testSpecialScalars()
    {
        BLS12_381G2Point g = BLS12_381G2.getGenerator();
        BigInteger[] scalars = { BigInteger.ZERO, BigInteger.ONE, BigInteger.valueOf(2), BigInteger.valueOf(3),
            R.subtract(BigInteger.ONE), R, R.add(BigInteger.ONE), R.add(BigInteger.valueOf(2)),
            BigInteger.ONE.shiftLeft(255), BigInteger.ONE.shiftLeft(256).subtract(BigInteger.ONE),
            BigInteger.valueOf(-1), R.negate().add(BigInteger.valueOf(5)) };
        for (int i = 0; i < scalars.length; ++i)
        {
            check(g, scalars[i], true);
        }
        assertTrue(g.multiply(R).isInfinity());
        assertTrue(g.constantTimeMultiply(R).isInfinity());
        assertEquals(g.doublePoint(), g.multiply(R.add(BigInteger.valueOf(2))));
    }

    public void testRandomScalarsInSubgroup()
    {
        BLS12_381G2Point g = BLS12_381G2.getGenerator().multiply(new BigInteger(255, RNG));
        for (int i = 0; i < 6; ++i)
        {
            check(g, new BigInteger(256, RNG), true);
        }
    }

    public void testPointsOutsideSubgroup()
    {
        for (int i = 0; i < 3; ++i)
        {
            BLS12_381G2Point p = nonSubgroupPoint();
            check(p, new BigInteger(256, RNG), true);
            check(p, R, true);
            // Cofactor clearing, the hash-to-curve use of multiply.
            check(p, BLS12_381G2.H_EFF, false);
        }
    }

    public void testInfinity()
    {
        assertTrue(BLS12_381G2Point.INFINITY.multiply(BigInteger.TEN).isInfinity());
        assertTrue(BLS12_381G2Point.INFINITY.constantTimeMultiply(BigInteger.TEN).isInfinity());
    }
}
