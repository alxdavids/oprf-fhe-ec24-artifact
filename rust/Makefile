CARGO=cargo
RELEASE=--release
MODE=fast
VENUE=ec
THREADS=64

.PHONY: build test run bench bench-paper

build:
	${CARGO} build ${RELEASE}

test:
	${CARGO} test ${RELEASE} -- --nocapture

run:
	${CARGO} run ${RELEASE}

bench:
	${CARGO} bench > benchmarks-${VENUE}-tMAX-${MODE}.txt

bench-paper:
	RAYON_NUM_THREADS=${THREADS} ${CARGO} bench > benchmarks-${VENUE}-t${THREADS}-${MODE}.txt ; \
	RAYON_NUM_THREADS=1 ${CARGO} bench > benchmarks-${VENUE}-t1-${MODE}.txt ; \
	
