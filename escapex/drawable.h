
#ifndef __DRAWABLE_H
#define __DRAWABLE_H

/* abstract interface to drawable things.
   think of these like GUI windows. */
struct Drawable {
  /* draw yourself on the surface.
     don't clear it, because someone
     else may have drawn before you. */
  virtual void draw() = 0;

  /* be notified that the screen size has changed. 
     screen->w and screen->h will have been updated
     with the new size. */
  virtual void screenresize() = 0;

  static void init();
};

extern Drawable *nodraw;

#endif

