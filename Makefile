

all: phenotype-code

phenotype-code:
	cd phenotype-schema && \
	protoc \
	-I . -I ../ga4gh-schemas/src/main/proto/ \
	--python_out=../ \
	phenotype.proto
