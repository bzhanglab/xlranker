<script lang="ts">
  import * as Select from "$lib/components/ui/select/index.js";
  import * as HoverCard from "$lib/components/ui/hover-card/index.js";
  import { Button, buttonVariants } from "$lib/components/ui/button/index.js";
  import FileUpload from "$lib/components/FileUpload.svelte";
  import { Loader } from "@lucide/svelte";
  import * as Sheet from "$lib/components/ui/sheet/index.js";
  let peptide_file: File[];
  let net_sep = "Tab";
</script>

<div class="px-4">
  <h3>Peptide Network</h3>
  Please upload a <Sheet.Root>
    <Sheet.Trigger
      class="decoration-dashed underline hover:decoration-solid font-semibold text-sm"
    >
      Peptide Network
    </Sheet.Trigger>
    <Sheet.Content>
      <Sheet.Header>
        <Sheet.Title>Peptide Network</Sheet.Title>
        <Sheet.Description>
          A text file containing peptide pairs that are observed to be
          cross-linked.
        </Sheet.Description>
      </Sheet.Header>

      <div class="grid flex-1 auto-rows-min gap-6 px-4 text-sm">
        <p>
          The peptide network is the list of peptides that were observed to be
          cross-linked. Each peptide should be the amino acid sequence (e.g. <code
            >KGMFR</code
          >).
        </p>
        <p>
          The file should have two columns, where each column represents a
          peptide. The order of the peptides does not matter and duplicated rows
          are removed automatically.
        </p>
        <p>
          After uploading, XLRanker will verify the format and alert of any
          potential issues.
        </p>

        <h5 class="font-semibold">Example</h5>
        <i>This example is tab-separated</i>
        <pre
          id="example"
          class="whitespace-pre-wrap bg-gray-100 p-4 rounded-md">
EQLAKLMATLR	IIGLDQVAGM
EQLAKLMATLR	KGMFR
EQLAKLMATLR	KGMFR
EQLDNQLDAK	SLNLKHIK
EQLDNQLDAY	SLNLKHIK
</pre>
        <h5 class="font-semibold">Supported Formats</h5>
        <ul class="list-disc list-inside">
          <li>Text files (.txt)</li>
          <li>Comma-separated values (.csv)</li>
          <li>Tab-separated values (.tsv)</li>
        </ul>

        <p>
          Specify the character separating by adjusting the <b
            >Column Separator</b
          > parameter.
        </p>
      </div>
      <Sheet.Footer>
        <Sheet.Close class={buttonVariants({ variant: "outline" })}
          >Close</Sheet.Close
        >
      </Sheet.Footer>
    </Sheet.Content>
  </Sheet.Root> file containing peptide pairs, one pair per line. Each line should
  contain two peptides separated by a tab or comma character.

  <div class="grid gap-2">
    <h3>Upload Network</h3>
    <FileUpload
      bind:files={peptide_file}
      accepted_formats={[".txt", ".csv", ".tsv"]}
    />
    <h4>Column Separator</h4>
    <Select.Root type="single" bind:value={net_sep}>
      <Select.Trigger class="w-fit min-w-[8rem]">{net_sep}</Select.Trigger>
      <Select.Content>
        <Select.Item value="Tab">Tab</Select.Item>
        <Select.Item value="Comma">Comma</Select.Item>
        <Select.Item value="Space">Space</Select.Item>
      </Select.Content>
    </Select.Root>
  </div>

  <!-- centered button for submit -->
  <div class="flex justify-center">
    <Button>
      <Loader style="animation: spin 2s linear infinite;" />
      Submit
    </Button>
  </div>
</div>
