package org.bouncycastle.pqc.crypto.faest;

import org.bouncycastle.util.Longs;

/**
 * Word-level carry-less multiplication shared by the FAEST binary fields.
 * <p>
 * Polynomials over GF(2) are held little-endian in 64-bit words: bit {@code i} of word {@code w} is the
 * coefficient of {@code x^(64w + i)}. A 64 x 64 bit product is formed with integer multiplications on
 * operands masked to every fourth bit (so the carries an integer multiply produces stay out of the bits that
 * are kept), which needs no table and no secret-dependent branch or index, keeping the arithmetic constant
 * time as the bit-serial loops it replaces were. The low half of the product comes from the operands, the
 * high half from their bit reversals. Multi-word products combine those word products by Karatsuba.
 */
final class BFMul
{
    private BFMul()
    {
    }

    /**
     * The low 64 bits of the carry-less product of {@code x} and {@code y}.
     * <p>
     * Each of the 16 integer products sums at most 16 one-bit terms per output position, so a carry out of a
     * 4-bit group needs 16 terms, which only output positions 60 and above collect; those carries land at
     * position 64 or higher, outside the result.
     */
    static long lo(long x, long y)
    {
        long x0 = x & 0x1111111111111111L;
        long x1 = x & 0x2222222222222222L;
        long x2 = x & 0x4444444444444444L;
        long x3 = x & 0x8888888888888888L;

        long y0 = y & 0x1111111111111111L;
        long y1 = y & 0x2222222222222222L;
        long y2 = y & 0x4444444444444444L;
        long y3 = y & 0x8888888888888888L;

        long z0 = (x0 * y0) ^ (x1 * y3) ^ (x2 * y2) ^ (x3 * y1);
        long z1 = (x0 * y1) ^ (x1 * y0) ^ (x2 * y3) ^ (x3 * y2);
        long z2 = (x0 * y2) ^ (x1 * y1) ^ (x2 * y0) ^ (x3 * y3);
        long z3 = (x0 * y3) ^ (x1 * y2) ^ (x2 * y1) ^ (x3 * y0);

        z0 &= 0x1111111111111111L;
        z1 &= 0x2222222222222222L;
        z2 &= 0x4444444444444444L;
        z3 &= 0x8888888888888888L;

        return z0 | z1 | z2 | z3;
    }

    /**
     * The high 64 bits of the carry-less product of {@code x} and {@code y}, given their bit reversals
     * {@code xr} and {@code yr}: the low half of the reversed product holds coefficients 63..126 in reverse.
     */
    static long hi(long xr, long yr)
    {
        return Longs.reverse(lo(xr, yr)) >>> 1;
    }

    /** {@code zz[zOff..zOff+2) ^= x * y}. */
    static void mul1Acc(long x, long y, long[] zz, int zOff)
    {
        zz[zOff]     ^= lo(x, y);
        zz[zOff + 1] ^= hi(Longs.reverse(x), Longs.reverse(y));
    }

    /** {@code zz[zOff..zOff+n+1) ^= x * y} for an {@code n}-word {@code x} and a 1-word {@code y}. */
    static void mulNx1Acc(long[] x, int xOff, int n, long y, long[] zz, int zOff)
    {
        long yr = Longs.reverse(y);
        for (int i = 0; i < n; ++i)
        {
            long xi = x[xOff + i];
            zz[zOff + i]     ^= lo(xi, y);
            zz[zOff + i + 1] ^= hi(Longs.reverse(xi), yr);
        }
    }

    /** {@code zz[zOff..zOff+4) ^= x * y} for 2-word operands, by Karatsuba (3 word products). */
    static void mul2Acc(long[] x, int xOff, long[] y, int yOff, long[] zz, int zOff)
    {
        long x0 = x[xOff], x1 = x[xOff + 1];
        long y0 = y[yOff], y1 = y[yOff + 1];
        long x0r = Longs.reverse(x0), x1r = Longs.reverse(x1);
        long y0r = Longs.reverse(y0), y1r = Longs.reverse(y1);

        long l0 = lo(x0, y0), h0 = hi(x0r, y0r);
        long l2 = lo(x1, y1), h2 = hi(x1r, y1r);
        long l1 = lo(x0 ^ x1, y0 ^ y1) ^ l0 ^ l2;
        long h1 = hi(x0r ^ x1r, y0r ^ y1r) ^ h0 ^ h2;

        zz[zOff]     ^= l0;
        zz[zOff + 1] ^= h0 ^ l1;
        zz[zOff + 2] ^= h1 ^ l2;
        zz[zOff + 3] ^= h2;
    }

    /**
     * {@code zz[zOff..zOff+6) ^= x * y} for 3-word operands, by one-level Karatsuba over all word pairs
     * (6 word products).
     */
    static void mul3Acc(long[] x, int xOff, long[] y, int yOff, long[] zz, int zOff)
    {
        long x0 = x[xOff], x1 = x[xOff + 1], x2 = x[xOff + 2];
        long y0 = y[yOff], y1 = y[yOff + 1], y2 = y[yOff + 2];
        long x0r = Longs.reverse(x0), x1r = Longs.reverse(x1), x2r = Longs.reverse(x2);
        long y0r = Longs.reverse(y0), y1r = Longs.reverse(y1), y2r = Longs.reverse(y2);

        long d0l = lo(x0, y0), d0h = hi(x0r, y0r);
        long d1l = lo(x1, y1), d1h = hi(x1r, y1r);
        long d2l = lo(x2, y2), d2h = hi(x2r, y2r);

        // t1 = x0 y1 + x1 y0, t2 = x0 y2 + x1 y1 + x2 y0, t3 = x1 y2 + x2 y1
        long t1l = lo(x0 ^ x1, y0 ^ y1) ^ d0l ^ d1l;
        long t1h = hi(x0r ^ x1r, y0r ^ y1r) ^ d0h ^ d1h;
        long t2l = lo(x0 ^ x2, y0 ^ y2) ^ d0l ^ d1l ^ d2l;
        long t2h = hi(x0r ^ x2r, y0r ^ y2r) ^ d0h ^ d1h ^ d2h;
        long t3l = lo(x1 ^ x2, y1 ^ y2) ^ d1l ^ d2l;
        long t3h = hi(x1r ^ x2r, y1r ^ y2r) ^ d1h ^ d2h;

        zz[zOff]     ^= d0l;
        zz[zOff + 1] ^= d0h ^ t1l;
        zz[zOff + 2] ^= t1h ^ t2l;
        zz[zOff + 3] ^= t2h ^ t3l;
        zz[zOff + 4] ^= t3h ^ d2l;
        zz[zOff + 5] ^= d2h;
    }

    /** {@code zz[zOff..zOff+8) ^= x * y} for 4-word operands, by two levels of Karatsuba (9 word products). */
    static void mul4Acc(long[] x, int xOff, long[] y, int yOff, long[] zz, int zOff)
    {
        long x0 = x[xOff], x1 = x[xOff + 1], x2 = x[xOff + 2], x3 = x[xOff + 3];
        long y0 = y[yOff], y1 = y[yOff + 1], y2 = y[yOff + 2], y3 = y[yOff + 3];
        long x0r = Longs.reverse(x0), x1r = Longs.reverse(x1), x2r = Longs.reverse(x2), x3r = Longs.reverse(x3);
        long y0r = Longs.reverse(y0), y1r = Longs.reverse(y1), y2r = Longs.reverse(y2), y3r = Longs.reverse(y3);

        // L = (x0 + x1 X)(y0 + y1 X)
        long a0l = lo(x0, y0), a0h = hi(x0r, y0r);
        long a2l = lo(x1, y1), a2h = hi(x1r, y1r);
        long a1l = lo(x0 ^ x1, y0 ^ y1) ^ a0l ^ a2l;
        long a1h = hi(x0r ^ x1r, y0r ^ y1r) ^ a0h ^ a2h;
        long l0 = a0l, l1 = a0h ^ a1l, l2 = a1h ^ a2l, l3 = a2h;

        // H = (x2 + x3 X)(y2 + y3 X)
        long b0l = lo(x2, y2), b0h = hi(x2r, y2r);
        long b2l = lo(x3, y3), b2h = hi(x3r, y3r);
        long b1l = lo(x2 ^ x3, y2 ^ y3) ^ b0l ^ b2l;
        long b1h = hi(x2r ^ x3r, y2r ^ y3r) ^ b0h ^ b2h;
        long h0 = b0l, h1 = b0h ^ b1l, h2 = b1h ^ b2l, h3 = b2h;

        // M = ((x0 + x2) + (x1 + x3) X)((y0 + y2) + (y1 + y3) X) - L - H
        long u0 = x0 ^ x2, u1 = x1 ^ x3, v0 = y0 ^ y2, v1 = y1 ^ y3;
        long u0r = x0r ^ x2r, u1r = x1r ^ x3r, v0r = y0r ^ y2r, v1r = y1r ^ y3r;
        long c0l = lo(u0, v0), c0h = hi(u0r, v0r);
        long c2l = lo(u1, v1), c2h = hi(u1r, v1r);
        long c1l = lo(u0 ^ u1, v0 ^ v1) ^ c0l ^ c2l;
        long c1h = hi(u0r ^ u1r, v0r ^ v1r) ^ c0h ^ c2h;
        long m0 = c0l ^ l0 ^ h0, m1 = c0h ^ c1l ^ l1 ^ h1, m2 = c1h ^ c2l ^ l2 ^ h2, m3 = c2h ^ l3 ^ h3;

        zz[zOff]     ^= l0;
        zz[zOff + 1] ^= l1;
        zz[zOff + 2] ^= l2 ^ m0;
        zz[zOff + 3] ^= l3 ^ m1;
        zz[zOff + 4] ^= h0 ^ m2;
        zz[zOff + 5] ^= h1 ^ m3;
        zz[zOff + 6] ^= h2;
        zz[zOff + 7] ^= h3;
    }

    /**
     * Reduce the {@code len}-word product {@code zz} modulo {@code x^(64n) + x^k1 + x^k2 + x^k3 + 1}
     * (with {@code 64 > k1 > k2 > k3 > 0} and {@code n >= 2}) into {@code dst[dOff..dOff+n)}. {@code zz} is
     * overwritten.
     */
    static void reduce(long[] zz, int len, int n, int k1, int k2, int k3, long[] dst, int dOff)
    {
        for (int i = len - 1; i >= n; --i)
        {
            long h = zz[i];
            zz[i - n]     ^= h ^ (h << k1) ^ (h << k2) ^ (h << k3);
            zz[i - n + 1] ^= (h >>> (64 - k1)) ^ (h >>> (64 - k2)) ^ (h >>> (64 - k3));
        }
        System.arraycopy(zz, 0, dst, dOff, n);
    }
}
