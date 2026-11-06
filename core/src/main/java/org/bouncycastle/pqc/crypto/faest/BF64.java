package org.bouncycastle.pqc.crypto.faest;

import org.bouncycastle.util.Longs;
import org.bouncycastle.util.Pack;

/**
 * GF(2^64) arithmetic for FAEST v2.0 universal hashing.
 * <p>
 * Reduction polynomial: x^64 + x^4 + x^3 + x + 1 (modulus low-word 0x1b).
 * Elements are bare {@code long} values (little-endian when serialized).
 * <p>
 * faest-ref source of truth: {@code fields.c} (bf64_*).
 */
final class BF64
{
    /** Reduction-polynomial low word. faest-ref: fields.c:15. */
    static final long MODULUS = (1L << 4) | (1L << 3) | (1L << 1) | 1L;

    static final int BYTES = 8;

    private BF64()
    {
    }

    /**
     * {@code a * b} in GF(2^64): one carry-less word product (see {@link BFMul}), then
     * the reduction. Constant time. faest-ref: {@code bf64_mul}.
     */
    static long mul(long a, long b)
    {
        long lo = BFMul.lo(a, b);
        long hi = BFMul.hi(Longs.reverse(a), Longs.reverse(b));
        // x^64 = x^4 + x^3 + x + 1; the bits shifted past x^63 by that fold are folded once more
        long s = (hi >>> 63) ^ (hi >>> 61) ^ (hi >>> 60);
        return lo ^ hi ^ (hi << 1) ^ (hi << 3) ^ (hi << 4) ^ s ^ (s << 1) ^ (s << 3) ^ (s << 4);
    }

    /** Load 8 little-endian bytes into a long. */
    static long load(byte[] src, int srcOff)
    {
        return Pack.littleEndianToLong(src, srcOff);
    }

    /** Store a long as 8 little-endian bytes. */
    static void store(byte[] dst, int dstOff, long src)
    {
        Pack.longToLittleEndian(src, dst, dstOff);
    }
}
