Ports from Airwindows.

# MODS

I'm removing some things both to make the code tighter and to ease porting.

## Floating point dithers

I'm commenting out all floating point to floating point dithers.

Rationale:

-IEEE Floating point doesn't truncate, but round to the nearest, and the noise
 level is scaled up or down in value with the exponent, so weak signals have
 the same truncation noise to signal ratio than louder signals. It's not so on
 integers.

-I null tested on console 6 and channel 9 On float (32 bit). Using JS bit meter
 there was noise starting at the 27th bit. Whatever the dither is doing I don't
 think it's extremely important.

-I'm going to use Reaper, which uses processDoubleReplacing. Handling truncation
 from 80 to 64 bit floats seems excessive.

This eases my porting.

In case I'm wrong. I left what I was doing disabled by the C preprocessor, see
e.g. channel9 or console6.

## Denormal routines

I also comment out the denormal prevention, as this is handled by Juce via
processor register settings.
