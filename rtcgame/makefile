all: rtcgame.js sim-test.js

rtcgame.js: domutil.ts netutil.ts net.ts sim.ts play.ts
	tsc --strict -t es2020 --out $@ $^

sim-test.js: sim.ts arcfour.ts ops-tester.ts sim-test.ts
	tsc --strict -t es2020 --out $@ $^

arcfour-test.js: arcfour.ts arcfour-test.ts
	tsc --strict -t es2020 --out $@ $^

controller-test.js: controller.ts controller-test.ts
	tsc --strict -t es2020 --out $@ $^
