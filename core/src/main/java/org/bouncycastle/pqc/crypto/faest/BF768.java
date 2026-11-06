package org.bouncycastle.pqc.crypto.faest;

import org.bouncycastle.util.Pack;

/**
 * GF(2^768) arithmetic for FAEST v2.0 universal-hashing accumulators.
 * <p>
 * Reduction polynomial: x^768 + x^19 + x^17 + x^4 + 1 (modulus low-word 0xA0011).
 * Elements are stored as twelve 64-bit limbs little-endian.
 * <p>
 * Asymmetric multiplication: BF768 &times; BF256 &rarr; BF768. Mirrors
 * {@code bf768_mul_256} in {@code faest-ref/fields.c:891}.
 * <p>
 * faest-ref source of truth: {@code fields.c} (bf768_*).
 */
final class BF768
{
    /** Reduction-polynomial low word. faest-ref: fields.c:27. */
    static final long MODULUS = (1L << 19) | (1L << 17) | (1L << 4) | 1L;

    static final int LIMBS = 12;
    static final int BYTES = 96;

    private BF768()
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
     * {@code dst := a * b} where {@code a} lives in GF(2^768) and {@code b} is a
     * GF(2^256) element read as a polynomial of degree below 256: each 256-bit chunk of
     * {@code a} is multiplied by {@code b} by Karatsuba (see {@link BFMul}), and the
     * product is reduced once against the BF768 modulus. Constant time. faest-ref:
     * {@code bf768_mul_256}.
     */
    static void mul256(long[] dst, int dOff,
                       long[] a, int aOff,
                       long[] b, int bOff)
    {
        long[] zz = new long[LIMBS + 4];
        for (int i = 0; i < LIMBS; i += 4)
        {
            BFMul.mul4Acc(a, aOff + i, b, bOff, zz, i);
        }
        BFMul.reduce(zz, LIMBS + 4, LIMBS, 19, 17, 4, dst, dOff);
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
