
// Generate PDF files.
//

#ifndef _CC_LIB_PDF_H
#define _CC_LIB_PDF_H

#include <stdint.h>
#include <stdio.h>
#include <optional>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <cstdint>

#include "image.h"
#include "hashing.h"

// Allows for quick generation of simple PDF documents.
// This is useful for producing easily printed output from C code, where
// advanced formatting is not required
//
// Note: All coordinates/sizes are in points (1/72 of an inch).
// All coordinates are based on 0,0 being the bottom left of the page.
// All colors are specified as a packed 32-bit value - see @ref PDF_RGB.
// Text strings are interpreted as UTF-8 encoded, but only a small subset of
// characters beyond 7-bit ascii are supported (see @ref pdf_add_text for
// details).

static constexpr inline float PDF_INCH_TO_POINT(float inch) {
  return inch * 72.0f;
}

static constexpr inline float PDF_MM_TO_POINT(float mm) {
  return mm * 72.0f / 25.4f;
}

struct PDF {
private:

  // Internal PDF objects. You can ignore this, but it has to be first.
  enum ObjType {
    OBJ_none,
    OBJ_info,
    OBJ_stream,
    OBJ_font,
    OBJ_page,
    OBJ_bookmark,
    OBJ_outline,
    OBJ_catalog,
    OBJ_pages,
    OBJ_image,
    OBJ_link,

    OBJ_count,
  };

  struct Object {
    Object(ObjType type) : type(type) {}
    virtual ~Object() {}

    ObjType type = OBJ_none;
    // PDF output index
    int index = 0;
    // Byte position within the output file
    int offset = 0;
    // Previous and next objects of this same type.
    Object *prev = nullptr;
    Object *next = nullptr;
  };

  class StreamObj;

public:

  // PDF has 14 built-in fonts.
  enum BuiltInFont {
    HELVETICA,
    HELVETICA_BOLD,
    HELVETICA_OBLIQUE,
    HELVETICA_BOLD_OBLIQUE,
    COURIER,
    COURIER_BOLD,
    COURIER_OBLIQUE,
    COURIER_BOLD_OBLIQUE,
    TIMES_ROMAN,
    TIMES_BOLD,
    TIMES_ITALIC,
    TIMES_BOLD_ITALIC,
    SYMBOL,
    ZAPF_DINGBATS,
  };

  // Metadata to be inserted into the header of the output PDF.
  // Because these fields must be null-terminated, the maximum
  // length is actually 63.
  struct Info {
    // Software used to create the PDF.
    char creator[64] = {};
    char producer[64] = {};
    // The title of the PDF (typically displayed in the
    // window bar when viewing).
    char title[64] = {};
    char author[64] = {};
    char subject[64] = {};
    // The date the PDF was created.
    char date[64] = {};
  };

  // Page is the only object that you need to think about.
  // Add them with AppendNewPage.
  struct Page : public Object {
    // Dimensions of the page, in points.
    float Width() const { return width; }
    float Height() const { return height; }

    void SetSize(float width, float height);

  private:
    // Created by AppendNewPage.
    Page() : Object(OBJ_page) {}
    float width = 0.0f;
    float height = 0.0f;
    std::vector<Object *> children;
    std::vector<Object *> annotations;
    friend class PDF;
  };

  // For fine-grained control over spacing.
  // Each substring has additional space after it (can be negative).
  // The final space does nothing and can take any value.
  //
  // Spacing is relative: It is simply scaled by the font size.
  //
  // Implementation note: Spacing here is nominally positive (larger
  // values mean more space) and in point units, but the underlying
  // TJ command uses -1/1000 points.
  using SpacedLine = std::vector<std::pair<std::string, float>>;

  // Just kidding. There's also FontObj!
  struct FontObj : public Object {
    FontObj() : Object(OBJ_font) {}

    // Returns the font's name as known to PDF.
    // For built-in fonts, this is something like "Helvetica-Bold".
    // For embedded fonts, it is something unique that need not
    // have anything to do with the font itself (e.g. "Font3").
    std::string BaseFont() const;

    // Return the kerning for the pair, if defined in the font. This
    // kerning value is the amount of additional space to add, but is
    // typically negative. Here, units are for the unscaled font
    // ("1pt"), so they are already appropriate for AddSpacedLine.
    //
    // Note that built-in fonts have no kerning tables, but you could
    // make your own.
    std::optional<double> GetKerning(int codepoint1, int codepoint2) const;

    // Kerns a single line of text.
    SpacedLine KernText(const std::string &text) const;

  private:
    // The font's index, like 3 for /F3.
    int font_index = 0;

    // For built-in fonts, the font.
    std::optional<BuiltInFont> builtin_font = std::nullopt;

    // For embedded fonts, the data stream.
    StreamObj *ttf = nullptr;
    // For embedded fonts, the table of widths.
    // These widths are scaled to "14pt".
    std::vector<uint16_t> widths;
    // Map from a pair of codepoints to the additional space.
    // This space is scaled to "1pt".
    std::unordered_map<std::pair<int, int>, double,
      Hashing<std::pair<int, int>>> kerning;

    const uint16_t *GetWidths() const;

    friend class PDF;
  };

  // A drawing operation within a path.
  //  See PDF reference for detailed usage.
  struct PathOp {
    char op;  /*!< Operation command. Possible operators are: m = move to, l =
                 line to, c = cubic bezier curve with two control points, v =
                 cubic bezier curve with one control point fixed at first
                 point, y = cubic bezier curve with one control point fixed
                 at second point, h = close path */
    float x1; /*!< X offset of the first point. Used with: m, l, c, v, y */
    float y1; /*!< Y offset of the first point. Used with: m, l, c, v, y */
    float x2; /*!< X offset of the second point. Used with: c, v, y */
    float y2; /*!< Y offset of the second point. Used with: c, v, y */
    float x3; /*!< X offset of the third point. Used with: c */
    float y3; /*!< Y offset of the third point. Used with: c */
  };

  // US Letter page.
  static constexpr float PDF_LETTER_WIDTH = PDF_INCH_TO_POINT(8.5f);
  static constexpr float PDF_LETTER_HEIGHT = PDF_INCH_TO_POINT(11.0f);

  static constexpr float PDF_A4_WIDTH = PDF_MM_TO_POINT(210.0f);
  static constexpr float PDF_A4_HEIGHT = PDF_MM_TO_POINT(297.0f);

  static constexpr float PDF_A3_WIDTH = PDF_MM_TO_POINT(297.0f);
  static constexpr float PDF_A3_HEIGHT = PDF_MM_TO_POINT(420.0f);

  // XXX to functions/constants. Maybe should switch to RGBA.

  // Pack R, G, B bytes into a 32-bit value.
  #define PDF_RGB(r, g, b)                                              \
    (uint32_t)((((r)&0xff) << 16) | (((g)&0xff) << 8) | (((b)&0xff)))

  // Pack A, R, G, B bytes into a 32-bit value in the order ARGB.
  #define PDF_ARGB(a, r, g, b)                                            \
    (uint32_t)(((uint32_t)((a)&0xff) << 24) | (((r)&0xff) << 16) |        \
               (((g)&0xff) << 8) | (((b)&0xff)))

  /*! Utility macro to provide bright red */
  #define PDF_RED PDF_RGB(0xff, 0, 0)

  /*! Utility macro to provide bright green */
  #define PDF_GREEN PDF_RGB(0, 0xff, 0)

  /*! Utility macro to provide bright blue */
  #define PDF_BLUE PDF_RGB(0, 0, 0xff)

  /*! Utility macro to provide black */
  #define PDF_BLACK PDF_RGB(0, 0, 0)

  /*! Utility macro to provide white */
  #define PDF_WHITE PDF_RGB(0xff, 0xff, 0xff)

  /*!
   * Utility macro to provide a transparent color
   * This is used in some places for 'fill' colors, where no fill is required
   */
  #define PDF_TRANSPARENT (uint32_t)(0xffu << 24)

  // Constructor. Give width and height of the page, and optional
  // header info.
  PDF(float width, float height, const std::optional<Info> &info);

  ~PDF();

  // If an operation fails, this gets the error code and a human-readable
  // error message.
  int GetErrCode() const;
  std::string GetErr() const;
  // Acknowledge an outstanding error.
  void ClearErr();

  // Set font to a font previously embedded with AddTTF.
  bool SetFont(const std::string &font_name);

  // Set font to one of the 14 built-ins.
  void SetFont(BuiltInFont font);

  // Must be a font previously added to this PDF object.
  void SetFont(FontObj *font);

  FontObj *GetCurrentFont() const { return current_font; }

  // Dimensions of the document, in points.
  float Width() const;
  float Height() const;

  // Returns null on failure.
  Page *AppendNewPage();

  // Retrieve a page by its number (starting from 1).
  // Note: The page must have already been created via AppendNewPage.
  // Returns null if the page is not found.
  Page *GetPage(int page_number);

  // Save the given pdf document to the supplied filename.
  // Pass the name of the file to store the PDF into (NULL for stdout)
  // Returns < 0 on failure, >= 0 on success
  // int pdf_save(struct pdf_doc *pdf, const char *filename);
  bool Save(const char *filename);

  // Add a line to the document. If page is null, then use the most
  // recently added page.
  void AddLine(float x1, float y1, float x2, float y2,
               float width, uint32_t color_rgb,
               Page *page = nullptr);

  void AddQuadraticBezier(
      // Start and end points.
      float x1, float y1, float x2, float y2,
      // Control point.
      float xq1, float yq1,
      float width,
      uint32_t color_rgb, Page *page = nullptr);

  void AddCubicBezier(
      // Start and end points.
      float x1, float y1, float x2, float y2,
      // Control points.
      float xq1, float yq1, float xq2, float yq2,
      float width,
      uint32_t color_rgb, Page *page = nullptr);

  void AddEllipse(
      // Center of the ellipse.
      float x, float y,
      // Radius on the x, y axes
      float xradius, float yradius,
      float width,
      uint32_t color, uint32_t fill_color,
      Page *page = nullptr);

  void AddCircle(float x, float y, float radius, float width,
                 uint32_t color,
                 uint32_t fill_color, Page *page = nullptr);

  void AddRectangle(float x, float y,
                    float width, float height, float border_width,
                    uint32_t color, Page *page = nullptr);

  void AddFilledRectangle(
    float x, float y, float width, float height,
    float border_width, uint32_t color_fill,
    uint32_t color_border, Page *page = nullptr);

  // Add a custom path as a series of ops (see PathOp).
  // Returns false if the ops are not understood.
  bool AddCustomPath(const std::vector<PathOp> &ops,
                     float stroke_width,
                     uint32_t stroke_color,
                     uint32_t fill_color,
                     Page *page = nullptr);

  // Returns false if the polygon is invalid (empty).
  bool AddPolygon(
      const std::vector<std::pair<float, float>> &points,
      float border_width, uint32_t color,
      Page *page = nullptr);

  // Returns false if the polygon is invalid (empty).
  // XXX I don't understand why this takes a single color
  // but has border_width?
  bool AddFilledPolygon(
    const std::vector<std::pair<float, float>> &points,
    float border_width, uint32_t color,
    Page *page = nullptr);

  // Barcodes.
  // https://en.wikipedia.org/wiki/Code_128
  bool AddBarcode128a(float x, float y, float width, float height,
                      const std::string &str, uint32_t color,
                      Page *page = nullptr);

  // https://en.wikipedia.org/wiki/Code_39
  bool AddBarcode39(float x, float y, float width, float height,
                    const std::string &str, uint32_t color,
                    Page *page = nullptr);

  // https://en.wikipedia.org/wiki/International_Article_Number
  bool AddBarcodeEAN13(float x, float y, float width, float height,
                       const std::string &str, uint32_t color,
                       Page *page = nullptr);

  // Encodes 12 digits.
  // https://en.wikipedia.org/wiki/Universal_Product_Code
  bool AddBarcodeUPCA(float x, float y, float width, float height,
                      const std::string &str, uint32_t color,
                      Page *page = nullptr);

  // Encodes 8 digits.
  // https://en.wikipedia.org/wiki/EAN-8
  bool AddBarcodeEAN8(float x, float y, float width, float height,
                      const std::string &str, uint32_t color,
                      Page *page = nullptr);

  // Encodes 12 digits with leading zeroes expected in certain
  // positions.
  // https://en.wikipedia.org/wiki/Universal_Product_Code
  bool AddBarcodeUPCE(float x, float y, float width, float height,
                      const std::string &str, uint32_t color,
                      Page *page = nullptr);

  bool AddText(const std::string &text,
               float size,
               // XXX baseline?
               float xoff, float yoff,
               uint32_t color, Page *page = nullptr);

  bool AddTextRotate(const std::string &text,
                     float size, float xoff, float yoff,
                     // In radians.
                     float angle,
                     uint32_t color,
                     Page *page = nullptr);

  // Text alignment for AddTextWrap.
  enum Alignment {
    PDF_ALIGN_LEFT,
    PDF_ALIGN_RIGHT,
    PDF_ALIGN_CENTER,
    // Align text in the center, with padding to fill the
    // available space.
    PDF_ALIGN_JUSTIFY,
    // Like PDF_ALIGN_JUSTIFY, except even short
    // lines will be fully justified
    PDF_ALIGN_JUSTIFY_ALL,
    // Fake alignment for only checking wrap height with
    // no writes
    PDF_ALIGN_NO_WRITE,
  };

  bool AddTextWrap(const std::string &text,
                   float size,
                   // XXX baseline?
                   float xoff, float yoff,
                   // in radians
                   float angle,
                   uint32_t color,
                   // Width of the box that text is written into.
                   float wrap_width,
                   Alignment alignment = PDF_ALIGN_LEFT,
                   // Height used, if non-null.
                   float *height = nullptr,
                   Page *page = nullptr);

  bool AddSpacedLine(const SpacedLine &line,
                     float size,
                     float xoff, float yoff,
                     uint32_t color,
                     float angle = 0.0f,
                     Page *page = nullptr);

  // Returns nullptr if the font has not yet been loaded
  // with SetFont (for built-in fonts) or AddTTF.
  FontObj *GetFontByName(const std::string &name) const;
  FontObj *GetBuiltInFont(BuiltInFont font) const;

  // Returns true upon success.
  bool GetTextWidth(const std::string &text,
                    float size,
                    // Out parameter.
                    float *text_width,
                    // If null, use current font.
                    FontObj *font = nullptr);

  // Add a bookmark to the document.
  // The page is the page to jump to (or nullptr for the most recent one).
  // The parent bookmark id, or -1 if this is a top-level bookmark.
  // Returns the non-negative bookmark id.
  int AddBookmark(const std::string &name, int parent, Page *page);

  bool AddLink(
      // The clickable rectangle.
      float x, float y,
      float width, float height,
      // Page that link should jump to.
      Page *target_page,
      // Point to place at the top left of the view.
      float target_x, float target_y,
      Page *page);

  enum class CompressionType {
    PNG,
    JPG_0,
    JPG_10,
    JPG_20,
    JPG_30,
    JPG_40,
    JPG_50,
    JPG_60,
    JPG_70,
    JPG_80,
    JPG_90,
    JPG_100,
  };

  bool AddImageRGB(float x, float y,
                   // If one of width or height is negative, then the
                   // value is determined from the other, preserving the
                   // aspect ratio.
                   float width, float height,
                   const ImageRGB &img,
                   CompressionType compression = CompressionType::PNG,
                   Page *page = nullptr);

  // TODO: Pretty easy to support 8-bit greyscale with ImageA.

  // TODO: Add support for RGBA images with actual alpha channels!

  // XXX In progress. TODO: return font name
  std::string AddTTF(const std::string &filename);

  static const char *ObjTypeName(ObjType t);

private:

  struct OutlineObj : public Object {
    OutlineObj() : Object(OBJ_outline) {}
  };

  struct BookmarkObj : public Object {
    BookmarkObj() : Object(OBJ_bookmark) {}
    Object *page = nullptr;
    std::string name;
    Object *parent = nullptr;
    std::vector<Object *> children;
  };

  struct StreamObj : public Object {
    StreamObj() : Object(OBJ_stream) {}
    Object *page = nullptr;
    std::string stream;
  };

  // Port note: This originally used the "stream" entry in the union.
  struct ImageObj : public Object {
    ImageObj() : Object(OBJ_image) {}
    Page *page = nullptr;
    std::string stream;
  };

  struct InfoObj : public Object {
    InfoObj() : Object(OBJ_info) {}
    Info info;
  };

  struct LinkObj : public Object {
    LinkObj() : Object(OBJ_link) {}
    // Page containing link.
    Page *page = nullptr;
    // Clickable rectangle.
    float llx = 0.0f;
    float lly = 0.0f;
    float urx = 0.0f;
    float ury = 0.0f;
    // Target page.
    Page *target_page = nullptr;
    // Target location.
    float target_x = 0.0f;
    float target_y = 0.0f;
  };

  struct NoneObj : public Object {
    NoneObj() : Object(OBJ_none) {}
  };

  struct CatalogObj : public Object {
    CatalogObj() : Object(OBJ_catalog) {}
  };

  struct PagesObj : public Object {
    PagesObj() : Object(OBJ_pages) {}
  };

  // Supported barcode guard patterns.
  enum GuardPattern {
    GUARD_NORMAL,
    GUARD_CENTRE,
    GUARD_SPECIAL,
    GUARD_ADDON,
    GUARD_ADDON_DELIN,
  };

  int SetErr(int errval, const char *buffer, ...);
  Object *pdf_get_object(int index);
  void pdf_append_object(Object *obj);
  // Perhaps can just be destructor, or static
  void pdf_object_destroy(Object *object);
  // Add a new-ly allocated object; takes ownership.
  Object *pdf_add_object_internal(Object *obj);
  // Convenience method that casts back to the derived class.
  template<class T> inline T *AddObject(T *t) {
    return (T*)pdf_add_object_internal(t);
  }
  void pdf_del_object(Object *obj);
  Object *pdf_find_first_object(int type);
  Object *pdf_find_last_object(int type);
  static int pdf_get_bookmark_count(const Object *obj);
  void pdf_add_stream(Page *page, std::string str);
  int pdf_save_file(FILE *fp);
  int pdf_save_object(FILE *fp, int index);

  bool pdf_text_point_width(const char *text,
                            ptrdiff_t text_len, float size,
                            const uint16_t *widths, float *point_width);

  bool pdf_add_text_spacing(const std::string &text, float size, float xoff,
                            float yoff, uint32_t color, float spacing,
                            float angle, Page *page);

  float pdf_barcode_128a_ch(float x, float y, float width, float height,
                            uint32_t color, int index, int code_len,
                            Page *page);
  bool pdf_barcode_39_ch(float x, float y, float char_width, float height,
                         uint32_t color, char ch, float *new_x, Page *page);

  bool pdf_barcode_eanupc_ch(float x, float y, float x_width,
                             float height, uint32_t color, char ch,
                             int set, float *new_x, Page *page);

  void pdf_barcode_eanupc_aux(float x, float y,
                              float x_width, float height,
                              uint32_t color, GuardPattern guard_type,
                              float *new_x, Page *page);

  bool pdf_add_image(ImageObj *image, float x, float y,
                     float width, float height, Page *page);

  ImageObj *pdf_add_raw_grayscale8(const uint8_t *data,
                                   uint32_t width,
                                   uint32_t height);

  bool pdf_add_png_data(float x, float y,
                        float display_width,
                        float display_height,
                        const uint8_t *png_data,
                        size_t png_data_length, Page *page);

  bool pdf_add_jpeg_data(float x, float y,
                         float display_width,
                         float display_height,
                         const uint8_t *jpeg_data,
                         size_t len,
                         Page *page);

  // XXX use std::string
  char errstr[128] = {};
  int errval = 0;
  std::vector<Object *> objects;

  float width = 0.0f;
  float height = 0.0f;

  FontObj *current_font = nullptr;
  int next_font_index = 1;

  // Indexed by font name. Font objects are not owned; they
  // are from the objects vector.
  std::unordered_map<std::string, FontObj *> embedded_fonts;
  std::unordered_map<BuiltInFont, FontObj *> builtin_fonts;

  Object *last_objects[OBJ_count] = {};
  Object *first_objects[OBJ_count] = {};
};

#endif

