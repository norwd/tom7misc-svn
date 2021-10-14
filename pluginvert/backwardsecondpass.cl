
// Second pass of error propagation from destination to the source.
// This acts only on the source layer and each node is independent,
// but we have to do this chunk-by-chunk because the derivative of
// the output depends on the chunk's transfer function.

// We also optionally support clipping the errors here.

// Expects the following defines:
//
// DERIVATIVE, the derivative of the transfer function for the chunk
//   given in terms of the function's output as usual.
// CHUNK_START, an integer giving the start of the chunk in the
//   SOURCE layer.
// CLIP_ERROR, a boolean. If true, error values are clamped.
// LARGE_ERROR, a positive float giving the maximum permissible post-
//   derivative error per node per round. This is not quite the same
//   as a limit on the learning rate, because it also prevents large
//   errors from being backpropagated. Ignored if CLIP_ERROR is false.
//   value I've used in the past: 1000000.0f

__kernel void BackwardSecondPass(
                  // Full src output, size layers[src].num_nodes.
                  __global const float *restrict src_output,
                  // Full src errors, parallel to src_output.
                  __global float *restrict src_error) {
  // index into chunk
  // PERF: We could avoid the constant offset by doing this
  // with global_work_offset?
  const int h_chunk = get_global_id(0);
  // index into src_output etc.
  const int h_global = CHUNK_START + h_chunk;
  const float out_h = src_output[h_global];
  const float weighted_error_sum = src_error[h_global];

  const float err = DERIVATIVE(out_h) * weighted_error_sum;

  src_error[h_global] =
  #if CLIP_ERROR
    fmax(-LARGE_ERROR, fmin(LARGE_ERROR, err))
  #else
    err
  #endif
    ;
}
