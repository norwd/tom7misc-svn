
#include <stdio.h>
#include <unistd.h>
#include <cstdint>
#include <vector>

#include "base/logging.h"
#include "base/stringprintf.h"
#include "pi/bcm2835.h"
#include "arcfour.h"
#include "timer.h"

using uint8 = uint8_t;
using uint32 = uint32_t;
using namespace std;

// Six bits at 22, 23, 24, 25, 26, 27
static constexpr int DEMUX_GPIO = 22;

static constexpr int GROUP_SEL_A = 17;
static constexpr int GROUP_SEL_B = 16;

// XXX make configurable
static constexpr bool USE_GROUP_A = true;

static void SetDemux(uint8 value) {
  CHECK_GE(value, 0) << value;
  CHECK_LT(value, 64) << value;

  {
    // This is hard-coded, but we need to set the group sel bits to
    // 10 or 01.
    bcm2835_gpio_fsel(GROUP_SEL_A, BCM2835_GPIO_FSEL_OUTP);
    bcm2835_gpio_set_pud(GROUP_SEL_A, BCM2835_GPIO_PUD_OFF);
    bcm2835_gpio_fsel(GROUP_SEL_B, BCM2835_GPIO_FSEL_OUTP);
    bcm2835_gpio_set_pud(GROUP_SEL_B, BCM2835_GPIO_PUD_OFF);

    constexpr uint32 GROUP_MASK =
      (1 << GROUP_SEL_A) | (1 << GROUP_SEL_B);
    if (USE_GROUP_A) {
      bcm2835_gpio_write_mask((1 << GROUP_SEL_A), GROUP_MASK);
    } else {
      bcm2835_gpio_write_mask((1 << GROUP_SEL_B), GROUP_MASK);
    }
  }
  
  for (int i = 0; i < 6; i++) {
    bcm2835_gpio_fsel(DEMUX_GPIO + i, BCM2835_GPIO_FSEL_OUTP);
    bcm2835_gpio_set_pud(DEMUX_GPIO + i, BCM2835_GPIO_PUD_OFF);
  }
    
  bcm2835_gpio_write_mask((value << DEMUX_GPIO), (0b111111 << DEMUX_GPIO));
  printf("Set demux to %d.\n", value);

  return;
}

int main(int argc, char **argv) {

  CHECK(bcm2835_init()) << "BCM Init failed!";

  CHECK(argc == 2) << "Usage: ./setdemux.exe value\n"
    "Run as root. value must be in [0, 7].\n";

  int value = atoi(argv[1]);
  printf("Setting demux to %d\n", value);
	 
  SetDemux(value);

  // blinkenlights demo
  #if 0
  uint8 v = 0;
  for (;;) {
    SetDemux(v & 0b11);
    usleep(50000);
    v++;
  }  
  #endif
  
  return 0;
}
