
// Second pass of the update weights phase.
//
// grad_sums   - sum of updates (gradients that tell us the desired
//               change in weight - these have the opposite sign
//               in some presentations) across all examples
// weights     - existing weights for the chunk; same size
// weights_aux - for ADAM, this is 2x the size of the above. for SGD, empty
//
// the loop is basically weights[i] += scale * grad_sums[i]
// where scale includes the learning rate and (1/num_examples), at least
// naively. In the case of ADAM, we do something more complicated and
// update weights_aux as well.
//
// weights here can actually be biases; we use the same kernel for each.

// Expects the following defines:
//
// EXAMPLES_PER_ROUND, the number of examples that contributed so that
//   we can compute the average.
// Either WEIGHT_UPDATE_SGD or WEIGHT_UPDATE_ADAM, which tells us what
//   method to use to compute the update (and how to interpret the aux
//   array).

// CLIPPING, if true, clips each grad to [-1,1] before applying any
//   update.
// CONSTRAIN, if true, ensures that the resulting weights are always
//   in [-constrain_max, constrain_max] (kernel parameter).

// TODO: Make configurable, although it seems that these are
// rarely tuned.
#define ADAM_B1 0.9f
#define ADAM_B2 0.999f
// This one apparently can use some tuning; I don't have any
// good intuitions though.
// XXXX
// #define ADAM_EPSILON 0.0000001f
#define ADAM_EPSILON 0.001f

__kernel void UpdateWeightsSecondPass(
                 // zero-based round number. PERF: if one-based, saves
                 // instruction
                 const int round_number,
                 const float learning_rate,
                 // Only used if CONSTRAIN is defined.
                 // PERF: Could be compile-time constant, but some
                 // complexity since it is different for weights and
                 // biases.
                 const float constrain_max,
                 // These two are the same size.
                 __global float *restrict grad_sums,
                 __global float *restrict chunk_weights,
                 // For SGD, empty. For ADAM, num_weights * 2
                 __global float *restrict chunk_weights_aux) {
  const int idx = get_global_id(0);

  // Average gradient over all examples.
  const float raw_grad = grad_sums[idx] * (1.0f / EXAMPLES_PER_ROUND);
  #if CLIPPING
    // fmin and fmax should reject nan, inf
    const float grad = fmax(-1.0f, fmin(1.0f, raw_grad));
  #else
    const float grad = raw_grad;
  #endif

    // compute the update u according to the method
  #if WEIGHT_UPDATE_SGD
    const float u = learning_rate * grad;
  #elif WEIGHT_UPDATE_ADAM
    // PERF: Skip _hat step when round is sufficiently large
    const int midx = idx * 2;
    const int vidx = idx * 2 + 1;
    const float m_prev = chunk_weights_aux[midx];
    const float v_prev = chunk_weights_aux[vidx];
    const float m_new = ADAM_B1 * m_prev + (1.0f - ADAM_B1) * grad;
    const float v_new = ADAM_B2 * v_prev + (1.0f - ADAM_B2) * (grad * grad);
    // TODO: Avoid nan poisoning here
    chunk_weights_aux[midx] = m_new;
    chunk_weights_aux[vidx] = v_new;
    const float m_hat = m_new / (1.0f - pow(ADAM_B1, round_number + 1));
    const float v_hat = v_new / (1.0f - pow(ADAM_B2, round_number + 1));
    const float u = learning_rate * (m_hat / (sqrt(v_hat) + ADAM_EPSILON));
  #else
    #error Weight update must be SGD or ADAM
  #endif

    // PERF -- generate a separate multiplier and value and use fma()
  #if CONSTRAIN
    chunk_weights[idx] = fmax(-constrain_max,
                              fmin(constrain_max,
                                   chunk_weights[idx] + u));
  #else
    chunk_weights[idx] += u;
  #endif
}
