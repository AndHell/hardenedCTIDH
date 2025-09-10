#include <string.h>
#include <assert.h>

#include "primes.h"
#include "random.h"

#include "../CTIDH/ctidh.h"

// WOMBat keygen
void ctidh_private(private_key *priv)
{
    memset(priv, 0, sizeof(private_key));

    uint8_t rnd;

    uint64_t batch_sumkeys = 0;
    for (uint32_t b = 0; b < primes_batches; b++)
    {
        random_wombats(priv->ells, batch_numkeys[b], batch_start[b], batch_stop[b], batch_sumkeys);
        batch_sumkeys += batch_numkeys[b];
    }

    for (uint32_t b = 0; b < WOMBATKEYS; b++)
    {
        randombytes(&rnd, 1);

        // 0 := Dummy (unused as we go dummy free!)
        // 1 := +
        // 2 := -
    
        priv->directions[b] = (rnd % 2) + 1;
    }
}
