# The primes, the zeros, and a machine that never goes negative

*The project explained the way you'd explain it at a blackboard to someone
who is quick, knows some physics, and doesn't speak the mathematical
dialect. Written 26 November 2026. Numbers and dead ends are in
`questions.md`. Nothing here is a proof; it is what we saw when we poked.*

## 1. Two lists that are each other's sound

Take the primes and mark them on a tape, but not at 2, 3, 5, 7. Mark
them at their logarithms, so the tape is a slide rule: equal distances
mean equal ratios. Put a click at each mark, and fainter clicks at the
powers (4, 8, 9...). A prime p clicks with loudness ln p / √p.

Now ask the question you'd ask of any recording: what tones is it made
of? The answer is a second list of numbers, the **zeros**: 14.13, 21.02,
25.01, 30.42... Those are the pitches. Play all those tones together and
they cancel everywhere except at the marks, where they pile up into the
clicks. Primes are the recording, zeros are its spectrum. Neither list
is more fundamental. Each is the other, heard differently.

## 2. What the famous guess says

A tone can be steady, or it can swell or fade as the tape runs. The
Riemann hypothesis says: **every tone is steady.** No swelling, no
fading, forever. (The "straight line" you hear about is just the
mathematician's way of drawing "steady".)

If even one tone swelled, the primes would have a lumpiness that grows
without limit. Everybody believes none does. Nobody can say why.

## 3. Weil's trick: turn it into an energy

Here's the clever part. Take any short pulse you like. Ask: how strongly
does this pulse ring at each of the tones? Square that and add over all
tones. If every tone is steady, each term is a square, so the total can
never be negative. If some tone swells, you can design a pulse that
makes the total negative.

So: **steady tones ⟺ this "energy" is never negative, for any pulse.**

That alone would be useless, because it needs the tones. But the same
number can be worked out entirely from the recording side:

> energy = (a smooth background term) − (the echoes)

The echoes are simple. For each prime: slide a copy of the pulse along
the tape by that prime's distance, see how much it overlaps the
original, multiply by the prime's loudness. That's all a prime does in
this machine. It asks, "how much do you look like yourself, moved by my
distance?"

The background is a smooth thing with no primes in it. It comes from
ordinary geometry of the number line, and it grows gently, like a
logarithm, with pitch.

So the puzzle is now: why can the echoes never beat the background?

## 4. The window makes it finite

A pulse of length L can only overlap a copy of itself if the copy is
moved by less than L. So a short pulse cannot feel the big primes at
all. A pulse as long as ln 14 feels 2, 3, 5, 7, 11, 13 and their small
powers, nothing else.

That is the window. It turns one infinite question into a ladder of
finite ones, each a matrix a laptop can handle. The guess is true if
every rung is fine. We can only ever check rungs.

## 5. Why there are levels that are almost exactly zero

Think of the machine as an instrument. It has natural pulses, each with
its own energy: its **levels**.

Now count. A pulse of length L, built from pitches up to some top pitch,
has only so many independent shapes: roughly length × bandwidth. That is
the same counting that tells you how many samples a sound of given
duration needs. And there are only so many tones below that top pitch.

Do the count for a window of size c, with top pitch 2πc, and you find
**about 2c more shapes than tones.** With that much freedom to spare you
can build pulses that ring at *none* of the tones. Silent at every one.
Their energy should be zero.

It isn't quite zero, because of the oldest rule in wave physics: a pulse
that is strictly limited in length cannot be strictly limited in pitch.
A little always leaks out past the top. That leak is tiny, and it gets
smaller ferociously fast as the window grows.

So the form has about 2c levels that are 0.000...001. Connes and Consani
found these and tied them to the classic theory of that leak. They are
the centre of the whole story: the machine comes absurdly close to zero,
from above, and never crosses.

## 6. Above the top pitch: hide in the gaps

Above 2πc the count flips. Now there are more tones than shapes, and a
pulse can't dodge them all. The best it can do is put its pitch in a
wide gap between two neighbouring tones and keep its skirts narrow.

That is exactly what the lowest ordinary levels do. Each one is a wave
parked in a wide gap between zeros. As the window grows it sinks slowly
in its gap, then a better gap takes over.

Connes's group measured the lowest ordinary level and noticed a rule of
thumb for it, which they could not explain. The gap picture says the
rule can't last: wide gaps keep turning up, and our numbers say the
level drops below their rule once the window passes about 250. That is
a prediction someone can check.

## 7. Now poke it: move a prime

So far the primes were in the right places. The interesting question is
what's special about the right places. So we moved the click for 7 to
7.007.

The machine goes negative at once. Move it the other way: negative
again. Plot how negative against where you put the click, and you get a
**tent**, with its sharp peak at exactly 7.

Why a tent and not a bowl? A ball in a bowl feels no force at the
bottom; the cost starts out flat. A tent is what you get when there are
*many* states sitting at the same level, here all those near-zero
pulses. A small disturbance acts on that whole family at once, and it
pushes some up and some down. Whichever way you move the prime, some
pulse in the family profits. Physicists know this from atoms: put an
atom with several equal levels in a weak field and they split, some up,
some down, in proportion to the field. Same thing. We counted: move one
prime and about half the near-zero levels go negative, every time.

And the steepness has an exact rule. A prime's hold on its position is
its loudness times how sharply the "overlap with your shifted self"
changes, for the worst pulse in the family. We checked that against the
machine to five digits. A cute by-product: for the *weight* of a prime
the same reasoning gives the numbers 1, 1.414, 1.618, 1.732..., which
are what you get for a particle hopping along a chain of 2, 3, 4, 5
sites. A prime acts on the near-zero pulses like a hop along a short
chain, as many sites as its distance fits into the window.

## 8. Move them all: the tape speed

Maybe no single prime can move, but the whole set could shift in some
coordinated way? We let the computer search every combination, blind.

Nothing works. But one combination came back ten times gentler than all
the rest: move every prime in proportion to its own distance along the
tape. That is simply **running the tape at a slightly different speed.**
All the clicks keep their pattern; the whole recording is stretched.

Why is that one nearly free? Because the echo part of the machine only
cares about the pattern. Stretch the tape and stretch your pulse with
it, and every echo is as it was. What does *not* go along is the
background. And here is the nice bit: the background grows like a
logarithm of pitch, and when you rescale the thing inside a logarithm,
the logarithm just shifts by a constant. So changing the tape speed
adds (or subtracts) nearly the same small amount to **every** level at
once.

That's what we saw. Speed the tape one way: all the near-zero levels
sink together. The other way: they all rise together, except one odd
pulse, which turned out to be a special shape the machine keeps pinned
at a fixed spot (mathematicians call it the pole; it is the shape that
carries "how many primes there are on average").

So: **the primes hold their own pattern. The background sets the tape
speed.**

## 9. What the tiny levels are

Now put sections 5 and 8 together.

Slow the tape a little and every near-zero level sinks by the same
amount. The machine goes negative the moment the smallest one reaches
zero. So the size of the smallest level is precisely how much the tape
speed could be off before the machine breaks.

We tested that directly on small windows: predicted and measured
breaking points agree to a tenth of a percent, and the exchange rate is
almost one for one.

So the mysterious tiny numbers have a plain meaning. **They are the
tolerance on the tape speed.** And since they are the leak from section
5, the tolerance shrinks about a hundred-thousand-fold for every step
the window grows. Using only the primes up to 8, the speed is already
pinned to thirty decimal places.

"Why is the machine never negative?" has become: "why is the tape
running at exactly the right speed, inside a tolerance that closes
without limit?"

## 10. The cousin that has been solved

There is a version of this whole story for curves instead of whole
numbers, and there the guess was proved long ago. We ran the same pokes.

Two differences jump out.

First, in the curve world every click sits on an evenly spaced grid.
All the distances are multiples of one step. For the real primes the
distances ln 2, ln 3, ln 5 share no common step at all. We met this
head on when we tried to rebuild the machine as a simple network of
sites and links: each prime wants its own chain, with its own step, and
no single set of sites serves them all. Small primes behave like a
pattern in pitch, primes near the window's size behave like a link
between its two ends, and nothing makes both simple at once.

Second, in the curve world the background is positive by itself, and
there are exactly free ways to slide the tones around. The solved case
has room, and safe ground to stand on. For the real primes the
background alone is *negative*, and every small change we tried breaks
the machine. The primes aren't sitting in a valley. They're balancing a
debt, exactly.

That may be the simplest way to see why the proof for the cousin has
never been carried across.

## 11. What we don't know

Everything above is a description of *how* the machine stays positive:
which parts hold what, how tightly, with what tolerance. It is not an
explanation of *why*.

The why would be an account of why the near-zero pulses exist at every
window: why the echoes of the primes cancel a negative background so
exactly that what is left is only the leak that wave physics forces.
The counting in section 5 says there is *room* for such pulses if the
tones are steady. It does not say the tones are steady. That step is
the Riemann hypothesis itself, and every road we have tried ends at
that same door.

What we can say is that the door looks a bit different from here. On
this side of it the question is not about infinitely many mysterious
tones. It is about one recording, one smooth background, and whether
the tape speed is exactly right.
