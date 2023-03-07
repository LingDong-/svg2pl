cmd:
	gcc -Wall -O3 examples/example.c -o svg2pl
cmdt:
	gcc -Wall -O3 examples/hershey_text.c -o svg2pl
testurl: cmd
	curl -s $(FILE) | ./svg2pl > out.svg
test1: cmd
	./svg2pl $(FILE) > $(FILE).o.svg
dltest:
	mkdir -p tests;\
	while read -r line; do \
		curl -s https://dev.w3.org/SVG/tools/svgweb/samples/svg-files/$$line > tests/$$line; \
	done < examples/urls.txt;
testfold: cmd
	while read -r line; do \
		./svg2pl tests/$$line > tests/$$line.o.svg; \
	done < examples/urls.txt;