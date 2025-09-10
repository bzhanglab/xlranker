<script lang="ts">
  import { Label } from "$lib/components/ui/label";

  let {
    files = $bindable(),
    multiple_upload = false,
    accepted_formats = [],
  } = $props();

  let file_text = $state("Click or drag file to this area to upload");
  let show_formats = $state(true);

  function updateFileText() {
    show_formats = false;
    file_text = files.length
      ? files.map((f) => f.name).join(", ")
      : "Click or drag file to this area to upload";
  }
  function handleFileChange(event: Event) {
    event.preventDefault();
    const target = event.target as HTMLInputElement;
    if (target.files) {
      files = Array.from(target.files);
    }

    updateFileText();
  }

  function handleDrop(event: DragEvent) {
    event.preventDefault();
    if (event.dataTransfer?.files) {
      files = Array.from(event.dataTransfer.files);
    }

    updateFileText();
  }

  function handleDragOver(event: DragEvent) {
    event.preventDefault();
  }
</script>

<!-- Card container for visual structure -->
<div class="flex justify-center">
  <div class="w-full max-w-md">
    <!-- Drag-and-drop upload area -->
    <div
      class="w-full content-center mx-auto border-1 text-gray-600 border-dashed border-gray-300 rounded-md p-6 text-center cursor-pointer hover:border-gray-500 hover:text-black transition"
      ondrop={handleDrop}
      ondragover={handleDragOver}
      role="button"
      tabindex="0"
    >
      <Label for="file-upload" class="cursor-pointer flex flex-col">
        <p class="text-md text-center">
          {file_text}{#if accepted_formats.length > 0 && show_formats}<br /><br
            /><i>Supported formats: {accepted_formats.join(", ")}</i>{/if}
        </p>
      </Label>
      <input
        id="file-upload"
        type="file"
        class="hidden"
        accept={accepted_formats.join(",")}
        onchange={handleFileChange}
        multiple={multiple_upload}
      />
    </div>
  </div>
</div>
