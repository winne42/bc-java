package org.bouncycastle.crypto.bls;

import java.math.BigInteger;

/**
 * Affine point on the BLS12-381 G2 curve {@code E: y^2 = x^3 + 4*(1+I)} over
 * Fp^2.
 * <p>
 * Points are held in affine coordinates. {@link #add} and {@link #doublePoint}
 * work affinely, one {@link Fp2Element#inverse()} each; the two scalar
 * multiplications run in Jacobian coordinates ({@code x = X/Z^2},
 * {@code y = Y/Z^3}) with mixed additions of the affine base point, so they
 * need a single inversion, to return to affine form at the end.
 * <p>
 * Two scalar-multiplication APIs are provided: {@link #multiply} is the
 * variable-time double-and-add — fast for public scalars (cofactor
 * clearing, subgroup checks) but unsafe for secret scalars; and
 * {@link #constantTimeMultiply} is a fixed-iteration ladder for use with
 * secret scalars (e.g. {@code sk * H(msg)} in BLS sign).
 */
public final class BLS12_381G2Point
{
    /** B coefficient: 4 * (1 + I). */
    public static final Fp2Element B = Fp2Element.of(4, 4);

    /** The point at infinity (identity element). */
    public static final BLS12_381G2Point INFINITY = new BLS12_381G2Point(null, null, true);

    private final Fp2Element x;
    private final Fp2Element y;
    private final boolean infinity;

    private BLS12_381G2Point(Fp2Element x, Fp2Element y, boolean infinity)
    {
        this.x = x;
        this.y = y;
        this.infinity = infinity;
    }

    /**
     * Constructs a G2 point from affine coordinates and verifies that
     * {@code (x, y)} satisfies the curve equation {@code y^2 = x^3 + 4*(1+I)}.
     */
    public static BLS12_381G2Point of(Fp2Element x, Fp2Element y)
    {
        Fp2Element lhs = y.square();
        Fp2Element rhs = x.square().mul(x).add(B);
        if (!lhs.equals(rhs))
        {
            throw new IllegalArgumentException("point not on BLS12-381 G2 curve E");
        }
        return new BLS12_381G2Point(x, y, false);
    }

    /**
     * Skips the curve-equation check. Used by callers that have already
     * verified the point is on the curve (e.g. immediately after iso_3
     * evaluation in hash-to-curve).
     */
    static BLS12_381G2Point ofUnchecked(Fp2Element x, Fp2Element y)
    {
        return new BLS12_381G2Point(x, y, false);
    }

    public boolean isInfinity()
    {
        return infinity;
    }

    public Fp2Element x()
    {
        return x;
    }

    public Fp2Element y()
    {
        return y;
    }

    public BLS12_381G2Point negate()
    {
        return infinity ? this : new BLS12_381G2Point(x, y.neg(), false);
    }

    public BLS12_381G2Point add(BLS12_381G2Point other)
    {
        if (infinity)
        {
            return other;
        }
        if (other.infinity)
        {
            return this;
        }
        if (x.equals(other.x))
        {
            if (y.equals(other.y))
            {
                return doublePoint();
            }
            // y == -other.y => P + (-P) = O
            return INFINITY;
        }
        // slope = (y2 - y1) / (x2 - x1)
        Fp2Element slope = other.y.sub(y).mul(other.x.sub(x).inverse());
        Fp2Element x3 = slope.square().sub(x).sub(other.x);
        Fp2Element y3 = slope.mul(x.sub(x3)).sub(y);
        return new BLS12_381G2Point(x3, y3, false);
    }

    public BLS12_381G2Point doublePoint()
    {
        if (infinity || y.isZero())
        {
            return INFINITY;
        }
        // slope = (3*x^2) / (2*y)  (curve has A = 0)
        Fp2Element threeXSquared = x.square().mulFp(BigInteger.valueOf(3));
        Fp2Element twoY = y.mulFp(BigInteger.valueOf(2));
        Fp2Element slope = threeXSquared.mul(twoY.inverse());
        Fp2Element x3 = slope.square().sub(x.mulFp(BigInteger.valueOf(2)));
        Fp2Element y3 = slope.mul(x.sub(x3)).sub(y);
        return new BLS12_381G2Point(x3, y3, false);
    }

    /**
     * Constant-time scalar multiplication, suitable for secret scalars
     * (e.g. {@code sk * H(msg)} in BLS sign).
     * <p>
     * Uses a fixed-iteration "double, conditionally add" ladder over 256
     * bits, with the conditional-add implemented as an array-indexed
     * select rather than an {@code if}. Both branches of every iteration
     * compute the same set of point operations regardless of the
     * scalar-bit value, so the per-bit timing does not depend on
     * the scalar.
     * <p>
     * <b>Caveats.</b> "Constant-time" here means the
     * scalar-bit-pattern-independent at the scalar-mult loop level.
     * The underlying Jacobian point ops still have data-dependent branches
     * for the infinity / equal-x cases of the mixed addition (the accumulator
     * is at infinity until the scalar's first set bit; the equal-x cases are
     * negligibly probable for random secret scalars on a prime-order
     * subgroup), the {@link BigInteger} field arithmetic is not constant
     * time, and the JVM
     * itself may introduce cache / GC / JIT timing variance that pure
     * Java cannot fully eliminate. Sufficient against a remote network
     * timing attacker on a typical workload; not a substitute for a
     * constant-time native implementation against a co-located
     * adversary with cache-line resolution.
     */
    public BLS12_381G2Point constantTimeMultiply(BigInteger scalar)
    {
        if (scalar == null)
        {
            throw new IllegalArgumentException("scalar must not be null");
        }
        if (scalar.signum() < 0)
        {
            return negate().constantTimeMultiply(scalar.negate());
        }
        if (infinity)
        {
            return INFINITY;
        }

        // Fixed-width scalar handling: read 256 bits regardless of the
        // actual bit length of scalar, so the iteration count carries no
        // information about the secret.
        final int FIXED_BITS = 256;
        Fp2Element[] r = JACOBIAN_INFINITY;
        Fp2Element[][] options = new Fp2Element[2][];
        for (int i = FIXED_BITS - 1; i >= 0; --i)
        {
            r = jacobianDouble(r);
            Fp2Element[] candidate = jacobianAddAffine(r, x, y);
            options[0] = r;
            options[1] = candidate;
            int bit = scalar.testBit(i) ? 1 : 0;
            r = options[bit];
        }
        return fromJacobian(r);
    }

    /**
     * Variable-time double-and-add scalar multiplication. Suitable for
     * non-secret scalars (e.g. cofactor clearing); not safe for secret
     * scalar use — see {@link #constantTimeMultiply}.
     */
    public BLS12_381G2Point multiply(BigInteger scalar)
    {
        if (scalar.signum() == 0 || infinity)
        {
            return INFINITY;
        }
        if (scalar.signum() < 0)
        {
            return negate().multiply(scalar.negate());
        }
        Fp2Element[] r = JACOBIAN_INFINITY;
        for (int i = scalar.bitLength() - 1; i >= 0; --i)
        {
            r = jacobianDouble(r);
            if (scalar.testBit(i))
            {
                r = jacobianAddAffine(r, x, y);
            }
        }
        return fromJacobian(r);
    }

    /** The point at infinity in Jacobian coordinates (any point with Z = 0). */
    private static final Fp2Element[] JACOBIAN_INFINITY = { Fp2Element.ONE, Fp2Element.ONE, Fp2Element.ZERO };

    /**
     * Jacobian doubling for a = 0 (dbl-2009-l, Explicit-Formulas Database):
     * 2M + 5S. Z = 0 (infinity) and Y = 0 (order 2) both give Z3 = 0, so
     * neither needs a branch.
     */
    private static Fp2Element[] jacobianDouble(Fp2Element[] p)
    {
        Fp2Element x1 = p[0], y1 = p[1], z1 = p[2];
        Fp2Element a = x1.square();
        Fp2Element b = y1.square();
        Fp2Element c = b.square();
        Fp2Element t = x1.add(b).square().sub(a).sub(c);
        Fp2Element d = t.add(t);
        Fp2Element e = a.add(a).add(a);
        Fp2Element f = e.square();
        Fp2Element x3 = f.sub(d.add(d));
        Fp2Element c8 = c.add(c);
        c8 = c8.add(c8);
        c8 = c8.add(c8);
        Fp2Element y3 = e.mul(d.sub(x3)).sub(c8);
        Fp2Element yz = y1.mul(z1);
        Fp2Element z3 = yz.add(yz);
        return new Fp2Element[]{ x3, y3, z3 };
    }

    /**
     * Mixed addition of a Jacobian point and an affine point (madd-2007-bl,
     * Explicit-Formulas Database): 7M + 4S. Handles {@code p} at infinity
     * and the doubling and inverse cases, which the formula does not.
     */
    private static Fp2Element[] jacobianAddAffine(Fp2Element[] p, Fp2Element x2, Fp2Element y2)
    {
        Fp2Element x1 = p[0], y1 = p[1], z1 = p[2];
        if (z1.isZero())
        {
            return new Fp2Element[]{ x2, y2, Fp2Element.ONE };
        }
        Fp2Element z1z1 = z1.square();
        Fp2Element u2 = x2.mul(z1z1);
        Fp2Element s2 = y2.mul(z1).mul(z1z1);
        Fp2Element h = u2.sub(x1);
        Fp2Element sy = s2.sub(y1);
        if (h.isZero())
        {
            // Same x: P + P is a doubling, P + (-P) is infinity.
            return sy.isZero() ? jacobianDouble(p) : JACOBIAN_INFINITY;
        }
        Fp2Element hh = h.square();
        Fp2Element i = hh.add(hh);
        i = i.add(i);
        Fp2Element j = h.mul(i);
        Fp2Element r = sy.add(sy);
        Fp2Element v = x1.mul(i);
        Fp2Element x3 = r.square().sub(j).sub(v.add(v));
        Fp2Element y1j = y1.mul(j);
        Fp2Element y3 = r.mul(v.sub(x3)).sub(y1j.add(y1j));
        Fp2Element z3 = z1.add(h).square().sub(z1z1).sub(hh);
        return new Fp2Element[]{ x3, y3, z3 };
    }

    /** Converts Jacobian {@code (X, Y, Z)} to the affine point {@code (X/Z^2, Y/Z^3)}. */
    private static BLS12_381G2Point fromJacobian(Fp2Element[] p)
    {
        if (p[2].isZero())
        {
            return INFINITY;
        }
        Fp2Element zInv = p[2].inverse();
        Fp2Element zInv2 = zInv.square();
        return new BLS12_381G2Point(p[0].mul(zInv2), p[1].mul(zInv2).mul(zInv), false);
    }

    public boolean equals(Object other)
    {
        if (!(other instanceof BLS12_381G2Point))
        {
            return false;
        }
        BLS12_381G2Point o = (BLS12_381G2Point)other;
        if (infinity || o.infinity)
        {
            return infinity == o.infinity;
        }
        return x.equals(o.x) && y.equals(o.y);
    }

    public int hashCode()
    {
        if (infinity)
        {
            return 0;
        }
        return x.hashCode() * 31 + y.hashCode();
    }

    public String toString()
    {
        if (infinity)
        {
            return "G2(infinity)";
        }
        return "G2(" + x + ", " + y + ")";
    }
}
