package org.bouncycastle.pqc.crypto.faest;

import org.bouncycastle.util.Pack;

/**
 * GF(2^384) arithmetic for FAEST v2.0 universal-hashing accumulators.
 * <p>
 * Reduction polynomial: x^384 + x^12 + x^3 + x^2 + 1 (modulus low-word 0x100D).
 * Elements are stored as six 64-bit limbs little-endian.
 * <p>
 * The only multiplication exposed is <em>asymmetric</em>: a BF384 element times
 * a BF128 element, returning a BF384 element. This mirrors {@code bf384_mul_128}
 * in {@code faest-ref/fields.c:727} &mdash; universal hashing never multiplies
 * two BF384 elements together.
 * <p>
 * faest-ref source of truth: {@code fields.c} (bf384_*).
 */
final class BF384
{
    /** Reduction-polynomial low word. faest-ref: fields.c:23. */
    static final long MODULUS = (1L << 12) | (1L << 3) | (1L << 2) | 1L;

    static final int LIMBS = 6;
    static final int BYTES = 48;

    private BF384()
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
     * {@code dst := a * b} where {@code a} lives in GF(2^384) and {@code b} is a
     * GF(2^128) element read as a polynomial of degree below 128: each 128-bit chunk of
     * {@code a} is multiplied by {@code b} by Karatsuba (see {@link BFMul}), and the
     * product is reduced once against the BF384 modulus. Constant time. faest-ref:
     * {@code bf384_mul_128}.
     */
    static void mul128(long[] dst, int dOff,
                       long[] a, int aOff,
                       long[] b, int bOff)
    {
        long[] zz = new long[LIMBS + 2];
        for (int i = 0; i < LIMBS; i += 2)
        {
            BFMul.mul2Acc(a, aOff + i, b, bOff, zz, i);
        }
        BFMul.reduce(zz, LIMBS + 2, LIMBS, 12, 3, 2, dst, dOff);
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
