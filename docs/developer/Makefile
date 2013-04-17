SHELL=/bin/bash
.PHONY: all tek__all clean cleancache distclean
.SUFFIXES:
all: tek__all

.tek_cache/developer.tex: developer.tex
	@echo -e "TEXPP\tdeveloper.tex"
	@mkdir -p ".tek_cache" >& /dev/null || true
	@texpp -i "developer.tex" -o ".tek_cache/developer.tex"

developer.pdf: .tek_cache/developer.tex .tek_cache/images/pinspec.png.pdf .tek_cache/images/high-level.png.pdf .tek_cache/images/classRegion.png.pdf .tek_cache/images/classSurface.png.pdf .tek_cache/images/classTally.png.pdf
	@echo -e "LATEX\tdeveloper.tex"
	@cd ".tek_cache" ; pdflatex -interaction=batchmode "developer.tex" >& /dev/null || pdflatex "developer.tex"
	@cd ".tek_cache" ; pdflatex -interaction=batchmode "developer.tex" >& /dev/null || pdflatex "developer.tex"
	@cp ".tek_cache/developer.pdf" "developer.pdf"

.tek_cache/images/classTally.png.pdf: .tek_cache/images/classTally.png
	@echo -e "CONVERT\timages/classTally.png"
	@mkdir -p ".tek_cache/images" >& /dev/null || true
	@convert "images/classTally.png" ".tek_cache/images/classTally.png.pdf"

.tek_cache/images/classTally.png: images/classTally.png
	@echo -e "IMGCP\timages/classTally.png"
	@mkdir -p ".tek_cache/images" >& /dev/null || true
	@cp "images/classTally.png" ".tek_cache/images/classTally.png"

.tek_cache/images/classSurface.png.pdf: .tek_cache/images/classSurface.png
	@echo -e "CONVERT\timages/classSurface.png"
	@mkdir -p ".tek_cache/images" >& /dev/null || true
	@convert "images/classSurface.png" ".tek_cache/images/classSurface.png.pdf"

.tek_cache/images/classSurface.png: images/classSurface.png
	@echo -e "IMGCP\timages/classSurface.png"
	@mkdir -p ".tek_cache/images" >& /dev/null || true
	@cp "images/classSurface.png" ".tek_cache/images/classSurface.png"

.tek_cache/images/classRegion.png.pdf: .tek_cache/images/classRegion.png
	@echo -e "CONVERT\timages/classRegion.png"
	@mkdir -p ".tek_cache/images" >& /dev/null || true
	@convert "images/classRegion.png" ".tek_cache/images/classRegion.png.pdf"

.tek_cache/images/classRegion.png: images/classRegion.png
	@echo -e "IMGCP\timages/classRegion.png"
	@mkdir -p ".tek_cache/images" >& /dev/null || true
	@cp "images/classRegion.png" ".tek_cache/images/classRegion.png"

.tek_cache/images/high-level.png.pdf: .tek_cache/images/high-level.png
	@echo -e "CONVERT\timages/high-level.png"
	@mkdir -p ".tek_cache/images" >& /dev/null || true
	@convert "images/high-level.png" ".tek_cache/images/high-level.png.pdf"

.tek_cache/images/high-level.png: images/high-level.png
	@echo -e "IMGCP\timages/high-level.png"
	@mkdir -p ".tek_cache/images" >& /dev/null || true
	@cp "images/high-level.png" ".tek_cache/images/high-level.png"

.tek_cache/images/pinspec.png.pdf: .tek_cache/images/pinspec.png
	@echo -e "CONVERT\timages/pinspec.png"
	@mkdir -p ".tek_cache/images" >& /dev/null || true
	@convert "images/pinspec.png" ".tek_cache/images/pinspec.png.pdf"

.tek_cache/images/pinspec.png: images/pinspec.png
	@echo -e "IMGCP\timages/pinspec.png"
	@mkdir -p ".tek_cache/images" >& /dev/null || true
	@cp "images/pinspec.png" ".tek_cache/images/pinspec.png"

tek__all: developer.pdf

cleancache:
	@rm -rf ".tek_cache" >& /dev/null || true

clean: cleancache
	@rm "developer.pdf" >& /dev/null || true
	@rm "developer.html" >& /dev/null || true

distclean: clean
	@rm "Makefile" >& /dev/null || true
	@rm "images/Makefile" >& /dev/null || true

