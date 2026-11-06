package org.bouncycastle.pqc.crypto.faest;

import org.bouncycastle.util.Pack;

/**
 * GF(2^576) arithmetic for FAEST v2.0 universal-hashing accumulators.
 * <p>
 * Reduction polynomial: x^576 + x^13 + x^4 + x^3 + 1 (modulus low-word 0x2019).
 * Elements are stored as nine 64-bit limbs little-endian.
 * <p>
 * Asymmetric multiplication: BF576 &times; BF192 &rarr; BF576. Mirrors
 * {@code bf576_mul_192} in {@code faest-ref/fields.c:806}.
 * <p>
 * faest-ref source of truth: {@code fields.c} (bf576_*).
 */
final class BF576
{
    /** Reduction-polynomial low word. faest-ref: fields.c:25. */
    static final long MODULUS = (1L << 13) | (1L << 4) | (1L << 3) | 1L;

    static final int LIMBS = 9;
    static final int BYTES = 72;

    private BF576()
    {
    }

    static void zero(long[] dst, int off)
    {
        for (int i = 0; i < LIMBS; i++)
        {
            dst[off + i] = 0L;
        }
    }

    static void one(long[] dst, int off)
    {
        dst[off] = 1L;
        for (int i = 1; i < LIMBS; i++)
        {
            dst[off + i] = 0L;
        }
    }

    static boolean equals(long[] a, int aOff, long[] b, int bOff)
    {
        for (int i = 0; i < LIMBS; i++)
        {
            if (a[aOff + i] != b[bOff + i])
            {
                return false;
            }
        }
        return true;
    }

    static void add(long[] dst, int dOff, long[] a, int aOff, long[] b, int bOff)
    {
        for (int i = 0; i < LIMBS; i++)
        {
            dst[dOff + i] = a[aOff + i] ^ b[bOff + i];
        }
    }

    static void addInPlace(long[] acc, int accOff, long[] x, int xOff)
    {
        for (int i = 0; i < LIMBS; i++)
        {
            acc[accOff + i] ^= x[xOff + i];
        }
    }

    /**
     * {@code dst := a * b} where {@code a} lives in GF(2^576) and {@code b} is a
     * GF(2^192) element read as a polynomial of degree below 192: each 192-bit chunk of
     * {@code a} is multiplied by {@code b} by Karatsuba (see {@link BFMul}), and the
     * product is reduced once against the BF576 modulus. Constant time. faest-ref:
     * {@code bf576_mul_192}.
     */
    static void mul192(long[] dst, int dOff,
                       long[] a, int aOff,
                       long[] b, int bOff)
    {
        long[] zz = new long[LIMBS + 3];
        for (int i = 0; i < LIMBS; i += 3)
        {
            BFMul.mul3Acc(a, aOff + i, b, bOff, zz, i);
        }
        BFMul.reduce(zz, LIMBS + 3, LIMBS, 13, 4, 3, dst, dOff);
    }

    static void load(long[] dst, int off, byte[] src, int srcOff)
    {
        for (int i = 0; i < LIMBS; i++)
        {
            dst[off + i] = Pack.littleEndianToLong(src, srcOff + i * 8);
        }
    }

    static void store(byte[] dst, int dstOff, long[] src, int off)
    {
        for (int i = 0; i < LIMBS; i++)
        {
            Pack.longToLittleEndian(src[off + i], dst, dstOff + i * 8);
        }
    }
}
