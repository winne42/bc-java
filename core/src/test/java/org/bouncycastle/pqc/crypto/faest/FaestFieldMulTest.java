package org.bouncycastle.pqc.crypto.faest;

import java.security.SecureRandom;

import org.bouncycastle.util.Arrays;
import org.bouncycastle.util.test.SimpleTest;

/**
 * Checks every FAEST field multiplication against a plain bit-by-bit reference: carry-less schoolbook
 * product, then reduction one bit at a time from the top. The operands include the all-ones words, which
 * maximise the carries the integer-multiply word product has to keep out of its result, single top bits, and
 * random values; outputs are also checked with {@code dst} aliasing {@code a}.
 */
public class FaestFieldMulTest
    extends SimpleTest
{
    private static final int ITERATIONS = 2000;

    private final SecureRandom random = new SecureRandom();

    public String getName()
    {
        return "FaestFieldMul";
    }

    public void performTest()
    {
        checkWordProduct();

        for (int it = 0; it < ITERATIONS; ++it)
        {
            checkSquare(128, new int[]{ 7, 2, 1 }, it);
            checkSquare(192, new int[]{ 7, 2, 1 }, it);
            checkSquare(256, new int[]{ 10, 5, 2 }, it);
            checkMixed(384, 128, new int[]{ 12, 3, 2 }, it);
            checkMixed(576, 192, new int[]{ 13, 4, 3 }, it);
            checkMixed(768, 256, new int[]{ 19, 17, 4 }, it);
            checkMul64(128, new int[]{ 7, 2, 1 }, it);
            checkMul64(192, new int[]{ 7, 2, 1 }, it);
            checkMul64(256, new int[]{ 10, 5, 2 }, it);

            long a = word(it), b = word(it + 3);
            long[] expected = reference(new long[]{ a }, new long[]{ b }, 64, new int[]{ 4, 3, 1 });
            isTrue("BF64.mul", BF64.mul(a, b) == expected[0]);
        }
    }

    private void checkWordProduct()
    {
        long[] edge = { 0L, 1L, -1L, Long.MIN_VALUE, 0x8888888888888888L, 0x1111111111111111L,
            0xFFFFFFFF00000000L, 0x00000000FFFFFFFFL };
        for (int i = 0; i < edge.length + ITERATIONS; ++i)
        {
            long x = i < edge.length ? edge[i] : random.nextLong();
            for (int j = 0; j < edge.length + 8; ++j)
            {
                long y = j < edge.length ? edge[j] : random.nextLong();
                long[] zz = clmul(new long[]{ x }, new long[]{ y });
                isTrue("BFMul.lo", BFMul.lo(x, y) == zz[0]);
                isTrue("BFMul.hi", BFMul.hi(Long.reverse(x), Long.reverse(y)) == zz[1]);
            }
        }
    }

    private void checkSquare(int bits, int[] taps, int it)
    {
        int n = bits / 64;
        long[] a = operand(n, it), b = operand(n, it + 1);
        long[] expected = reference(a, b, bits, taps);
        long[] dst = new long[n];
        switch (bits)
        {
        case 128: BF128.mul(dst, 0, a, 0, b, 0); break;
        case 192: BF192.mul(dst, 0, a, 0, b, 0); break;
        default:  BF256.mul(dst, 0, a, 0, b, 0); break;
        }
        isTrue("BF" + bits + ".mul", Arrays.areEqual(expected, dst));

        long[] alias = Arrays.clone(a);
        switch (bits)
        {
        case 128: BF128.mul(alias, 0, alias, 0, b, 0); break;
        case 192: BF192.mul(alias, 0, alias, 0, b, 0); break;
        default:  BF256.mul(alias, 0, alias, 0, b, 0); break;
        }
        isTrue("BF" + bits + ".mul aliased", Arrays.areEqual(expected, alias));
    }

    private void checkMixed(int bits, int smallBits, int[] taps, int it)
    {
        long[] a = operand(bits / 64, it), b = operand(smallBits / 64, it + 1);
        long[] expected = reference(a, b, bits, taps);
        long[] dst = new long[bits / 64];
        switch (bits)
        {
        case 384: BF384.mul128(dst, 0, a, 0, b, 0); break;
        case 576: BF576.mul192(dst, 0, a, 0, b, 0); break;
        default:  BF768.mul256(dst, 0, a, 0, b, 0); break;
        }
        isTrue("BF" + bits + ".mul" + smallBits, Arrays.areEqual(expected, dst));
    }

    private void checkMul64(int bits, int[] taps, int it)
    {
        int n = bits / 64;
        long[] a = operand(n, it);
        long b = word(it + 5);
        long[] expected = reference(a, new long[]{ b }, bits, taps);
        long[] dst = new long[n];
        switch (bits)
        {
        case 128: BF128.mul64(dst, 0, a, 0, b); break;
        case 192: BF192.mul64(dst, 0, a, 0, b); break;
        default:  BF256.mul64(dst, 0, a, 0, b); break;
        }
        isTrue("BF" + bits + ".mul64", Arrays.areEqual(expected, dst));
    }

    /** Edge-case words on a fixed rotation, random words otherwise. */
    private long word(int i)
    {
        switch (i % 5)
        {
        case 0:  return -1L;
        case 1:  return Long.MIN_VALUE;
        default: return random.nextLong();
        }
    }

    private long[] operand(int n, int it)
    {
        long[] v = new long[n];
        for (int i = 0; i < n; ++i)
        {
            v[i] = (it & 1) == 0 ? word(it + i) : random.nextLong();
        }
        return v;
    }

    private static long[] clmul(long[] x, long[] y)
    {
        long[] zz = new long[x.length + y.length];
        for (int i = 0; i < x.length * 64; ++i)
        {
            if (((x[i >>> 6] >>> (i & 63)) & 1L) != 0)
            {
                for (int j = 0; j < y.length * 64; ++j)
                {
                    if (((y[j >>> 6] >>> (j & 63)) & 1L) != 0)
                    {
                        int k = i + j;
                        zz[k >>> 6] ^= 1L << (k & 63);
                    }
                }
            }
        }
        return zz;
    }

    /** x * y reduced modulo x^bits + x^taps[0] + x^taps[1] + x^taps[2] + 1. */
    private static long[] reference(long[] x, long[] y, int bits, int[] taps)
    {
        long[] zz = clmul(x, y);
        for (int k = zz.length * 64 - 1; k >= bits; --k)
        {
            if (((zz[k >>> 6] >>> (k & 63)) & 1L) != 0)
            {
                zz[k >>> 6] ^= 1L << (k & 63);
                int base = k - bits;
                flip(zz, base);
                for (int t = 0; t < taps.length; ++t)
                {
                    flip(zz, base + taps[t]);
                }
            }
        }
        return Arrays.copyOfRange(zz, 0, bits / 64);
    }

    private static void flip(long[] zz, int k)
    {
        zz[k >>> 6] ^= 1L << (k & 63);
    }

    public static void main(String[] args)
    {
        runTest(new FaestFieldMulTest());
    }
}
