default : markers.pdf


# Note 23 Jan 2022 ... I think this conversion does not work
# (output scale is wrong). I think I probably printed SVGs manually
# from illustrator?
MAKEPDF=../../pactom/rsvg/rsvg-convert -f pdf
markers.pdf : page1.svg page2.svg page3.svg page4.svg page5.svg page6.svg page7.svg page8.svg page9.svg page10.svg page11.svg
	${MAKEPDF} -a --dpi-x=90 --dpi-y=90 --width=765 $^ -o $@
