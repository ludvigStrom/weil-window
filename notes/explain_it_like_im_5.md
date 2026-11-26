# Explain it like I'm 5

The whole project, from the start, in plain words. Written 26 November
2026. For the careful version with numbers, see `questions.md`,
`ladder.md`, `string.md`, `function_field.md` and `shape_brainstorm.md`.

A warning first. This is a hobby project. None of it is a proof of
anything. Some of it is certainly known to experts, some of it might be
new, and some of it might be wrong.

## 1. What is the puzzle?

Prime numbers are the numbers you cannot split: 2, 3, 5, 7, 11, 13...
They look scattered at random, but they are not random. There is a
hidden list of other numbers, called the **zeros**, that works like the
sheet music for the primes. If you know the zeros you can rebuild the
primes, and the other way round.

The most famous unsolved puzzle in mathematics, the **Riemann
hypothesis**, says: all those zeros stand in one perfectly straight line.

Nobody knows why that should be true. Everybody believes it.

## 2. The machine

A mathematician called André Weil found another way to say the same
thing. Imagine a machine. You feed it a wiggle (any small wave you like)
and it gives you back a number.

Weil showed: **the zeros are in a straight line exactly when the machine
never gives a negative number.** No wiggle can make it go below zero.

So the puzzle becomes: why can the machine never go negative? That is
the question this project plays with. We don't try to prove it. We poke
the machine to see how it behaves.

## 3. The window

The full machine is too big. So we only allow wiggles that fit inside a
**window**, like only looking at a street through a letterbox. A small
window only "sees" the small primes: a window of size 14 sees 2, 3, 5,
7, 11 and 13 and nothing else.

That makes the machine small enough for a computer. All our experiments
are done on windows.

## 4. The machine's notes

Like a guitar string, the machine has its own favourite wiggles, and
each one has a number: its **level**. If every level is zero or more,
the machine never goes negative.

What we see (and what Alain Connes and his colleagues found before us):

- Most levels are ordinary sizes, like 0.2 or 1.5.
- But a handful are **incredibly tiny**: 0.000...0001 with thirty or
  fifty zeros. Tiny, but never below zero.
- Each time the window grows by one, one more tiny level appears.

Those tiny levels are the heart of the mystery. The machine comes
*unbelievably* close to going negative, and never does.

## 5. Taking a prime away

First poke: what if we remove one prime, say 7, and run the machine
again?

It breaks. One level drops below zero. Each prime is doing a job: it
cancels a bit of negativity that would otherwise be there. We measured
how far the machine falls without each prime, and how soon after the
window passes the prime the trouble starts.

## 6. A picture to think with

We wanted a physical picture, something you could build in a lab in
your imagination. The one that fits best:

- **A gas in a box.** The window is the box. The tiny levels behave like
  the quietest states of a cold gas squeezed into it.
- **Mirrors.** Each prime acts like a faint mirror placed at a special
  spot. Waves in the box bounce off them and leave faint echoes. We
  found those echoes in the computer, where the picture said they should
  be.
- **A horizon.** At the edges of the box, things behave like the edge of
  a black hole or of an expanding universe. That part was already in the
  work of Connes and Moscovici.

It's only a picture. But it was good enough to suggest experiments, and
the experiments came out the way it predicted.

## 7. The ordinary levels live in the gaps

Second finding. Remember the zeros, the sheet music? They are not evenly
spaced. Sometimes there is a wider gap between two neighbours.

We found: **the lowest of the ordinary levels belongs to a wave that
hides in one of the wide gaps between zeros.** Like a child playing
hide and seek who picks the biggest gap between the trees.

As the window grows, the wave stays in one gap for a while, then jumps
to a better one.

Connes's group had measured this lowest ordinary level and noticed it
follows a simple rule of thumb, but said they couldn't explain it. Our
picture gives a reason for the wobbles around their rule ("which gap is
in use"). And it makes a prediction: for bigger windows than they
tried, their rule of thumb should stop working. That is something they
could check.

(One piece of this, a shortcut for computing the levels straight from
the zeros, turned out to be already published by someone else. We now
say so.)

## 8. Moving a prime: the tent

Third poke, and the best one. Instead of removing the prime 7, we
**nudged** it. Pretend 7 sits at 7.001, or 6.999.

The machine breaks immediately, in both directions. If you draw how far
it falls against where you put the 7, you get a **tent**: a sharp peak
exactly at 7, sloping down on both sides.

So each prime sits on the tip of a tent. Move it a hair either way and
the machine goes negative. Small primes like 2 and 3 are pinned the
hardest, and a bigger window pins every prime harder.

## 9. Moving all the primes at once: the zoom

Maybe one prime can't move, but could several move *together* in some
clever way and keep the machine happy?

We let the computer search every combination, without telling it what to
look for.

- No combination works. The true primes are the only arrangement nearby
  that the machine accepts.
- But one combination is **ten times gentler** than all the others, and
  the computer found it by itself: move every prime the same way, like
  stretching a rubber ruler evenly. A **zoom**.

That zoom is not a random thing. It is exactly the motion that Connes's
whole approach to the puzzle is built around. We didn't put it in; the
machine pointed at it.

## 10. Who holds what

Why is the zoom almost allowed? We took the machine apart to see.

The machine is made of two things: the primes, and a smooth
**background** that has nothing to do with primes.

- The **pattern** of the primes, where they sit relative to each other,
  is held in place by the primes themselves, very hard.
- The **overall size**, the unit of length, is held by the background.

Like beads on a rubber band lying on a table with a ruler printed on it.
The beads keep their spacing pattern by themselves. Only the printed
ruler says how stretched the band should be.

## 11. What the tiny levels are

And now the tiny levels from section 4 get a meaning.

If you zoom the primes *out* a little, all the tiny levels sink together.
The machine breaks when the smallest one reaches zero. We checked: the
size of the smallest tiny level is, almost exactly, **how much zoom the
primes can take before the machine breaks**.

So the tiny levels are not mysterious small numbers. They are the
**wiggle room in the primes' overall scale**. And that room shrinks
about 100,000 times for every step the window grows. With only the
primes up to 8, the scale is already fixed to 30 decimal places.

Said another way: the question "why does the machine never go negative?"
becomes "why are the primes at exactly the right size, with less and
less room to spare, forever?"

## 12. The cousins where the puzzle is solved

There is a cousin of this puzzle, about curves instead of numbers, where
the answer *is* known: Weil proved it in the 1940s. We ran the same
pokes there.

- There, the "primes" also sit on tents.
- But the arrangement has **neighbours**: you can slide things around in
  three exactly free ways and the machine stays happy. And the
  background there is friendly: it is positive all by itself.

For the real primes there are no free moves, and the background by
itself is *negative*. The primes have to cancel it exactly.

That might be a way to see why the proof for the cousin has never been
carried over. The cousin lives in a roomy place. The real primes are
balanced on a pin.

## 13. What we have not done

- We have not proved anything.
- We have not explained **why** the primes hold their pattern so
  rigidly. That is the real remaining mystery.
- Everything is computed for small windows (mostly up to 14 or 28, a few
  things up to 1000). The real puzzle is about windows growing forever.
- Every road we have tried ends at the same wall: "true for each window
  we can compute; the puzzle is whether it stays true forever."

## 14. The one-paragraph version

There is a machine that never gives a negative number if, and only if,
the most famous guess in mathematics is true. We poked it. Every prime
sits on the tip of a tent: nudge it and the machine breaks. The primes
hold their own pattern rigidly, and the only thing with any give is the
overall zoom, which is held by a smooth background instead. The
machine's famous almost-zero levels turn out to be exactly that give,
and it shrinks at a furious rate as you look at more primes. In the
cousin problem that has been solved, there is plenty of room; for the
real primes there is none. We don't know why.
